from sgevalviz.reader import Reader
from pandas import DataFrame, Series
import pandas as pd
import numpy as np

class FillDataHelper:
    def __init__(self, chromosomePath: str, groupType: str, reader: Reader):
        self.reader = reader
        self.chromosomePath = chromosomePath
        self.groupType = groupType

        self.path = reader.getDefinedChromosomePath(chromosomePath, groupType, "processedGtf")
        self.geneTranscriptPath = reader.getDefinedChromosomePath(chromosomePath, groupType, "geneTranscript")

        self.df: DataFrame = pd.read_csv(self.path)
        self.geneStringPath = reader.getDefinedChromosomePath(chromosomePath, groupType, "geneString")

        self.dfString = None

        self.dfGeneString = None
        self.geneStringCompletePath = reader.getDefinedChromosomeSingleGeneStringPath(chromosomePath)

        self.geneTranscriptDf: DataFrame = pd.read_csv(self.geneTranscriptPath)

        self.isForwardStrand = self.df.iloc[0]["is_forward_strand"]

        self.dfExon = None
        self.dfIntron = None
        self.dfNotIntronOrExon = None

    def sortBy(self, sortList: list[str]):
        self.df = self.df.sort_values(by=sortList).reset_index(drop=True)

    def getUniqueGeneString(self):
        return self.dfGeneString["gene_string"].unique()

    def getMasks(self, maskList: list[str]) -> list[Series]:
        return [self.df[mask] for mask in maskList]

    def getSubDfs(self, boolSeriesList: list[Series]) -> list[DataFrame]:
        return [self.df.loc[boolSeries] for boolSeries in boolSeriesList]

    def setIntronExonAndOthers(self):
        isIntron, isExon = self.getMasks(["is_intron", "is_exon"])
        isNotIntronOrExon = ~isIntron & ~isExon

        dfExon, dfIntron, dfNotIntronOrExon = self.getSubDfs([isExon, isIntron, isNotIntronOrExon])
        self.dfExon = dfExon
        self.dfIntron = dfIntron
        self.dfNotIntronOrExon = dfNotIntronOrExon

    def defineFirstLastExon(self):
        g = self.dfExon.groupby(["gene_id","transcript_id"])

        smaller = g.head(1).index
        bigger = g.tail(1).index

        firstExon = smaller if self.isForwardStrand else bigger
        lastExon = bigger if self.isForwardStrand else smaller

        self.dfExon.loc[firstExon, "is_first_exon"] = True
        self.dfExon.loc[lastExon, "is_last_exon"] = True

    def defineIntronRetentionExon(self):
        self.dfExon["exon_start_repeats"] = self.dfExon.duplicated("region_start", keep=False)
        self.dfExon["exon_end_repeats"] = self.dfExon.duplicated("region_end", keep=False)

        min_end = self.dfExon.groupby("region_start")["region_end"].transform("min")
        self.dfExon["non_smallest_exon_from_same_start_group"] = (self.dfExon["region_end"] != min_end)

        self.dfExon["is_intron_retention_exon"] = (self.dfExon["exon_start_repeats"] & self.dfExon["exon_end_repeats"] & self.dfExon["non_smallest_exon_from_same_start_group"])
        self.dfExon.drop(columns=["exon_start_repeats", "exon_end_repeats", "non_smallest_exon_from_same_start_group"], inplace=True)

    def dropLastIntron(self):
        dfLocal = self.dfIntron.groupby(['gene_id', 'transcript_id'])

        lastIntronId = dfLocal.tail(1).index
        self.dfIntron.drop(lastIntronId, inplace=True)

    def unifyDf(self):
        df = pd.concat([self.dfExon, self.dfIntron]).sort_index().reset_index(drop=True)
        df = df.sort_values(by=['gene_id', 'transcript_id', 'region_start']).reset_index(drop=True)

        df['next_exon_start'] = np.where(df['is_exon'], df['region_start'], np.nan)
        df['next_exon_start'] = df['next_exon_start'].bfill()
        df.loc[df['is_intron'], 'region_end'] = df.loc[df['is_intron'], 'next_exon_start'] - 1
        df.drop(columns='next_exon_start', inplace=True)

        self.df = pd.concat([df, self.dfNotIntronOrExon]).sort_index().reset_index(drop=True)
        self.df['region_end'] = pd.to_numeric(self.df['region_end'], downcast='integer', errors='coerce')

    def enrinchDf(self):
        self.setIntronExonAndOthers()
        self.defineFirstLastExon()
        self.defineIntronRetentionExon()
        self.dropLastIntron()
        self.unifyDf()

    def addCodons(self):
        codonNames = ["start_codon", "stop_codon"]

        for codonName in codonNames:
            codonDf = self.df.loc[self.df[f"is_{codonName}"]]
            codonDf = codonDf[['chromosome_identifier','gene_id','transcript_id','region_start']]
            codonDf.rename(columns={'region_start':f"{codonName}_init"},inplace=True)
            self.dfString = pd.merge(self.dfString, codonDf, on=["chromosome_identifier", "gene_id", "transcript_id"], how="left") 

        cols = ["start_codon_init", "stop_codon_init", "gene_string"]
        self.dfString[cols] = self.dfString[cols].fillna("")

        self.dfString["gene_string"] = (
            "|"
            + self.dfString["start_codon_init"].astype(str)
            + "|"
            + self.dfString["gene_string"].astype(str)
            + "|"
            + self.dfString["stop_codon_init"].astype(str)
            + "|"
        )

        self.dfString["predicted"] = False
        self.dfString["gene_predicted"] = False


    def generateGeneStringDf(self):
        dfString = self.df.loc[self.df['is_exon']].copy()
        dfString['gene_string'] = dfString['region_start'].astype(str) + ';' + dfString['region_end'].astype(str)

        self.dfString = dfString.groupby(['chromosome_identifier','gene_id', 'transcript_id']
        ).agg(
            min_pos=("region_start", "min"), 
            max_pos=("region_end", "max"), 
            gene_string=('gene_string', '/'.join), # join all gene_string values
            exon_qtty=('gene_string', 'size'), # count how many were merged
            intron_retention_qtty=('is_intron_retention_exon','sum') #count how many of the exons are intron retention
        ).reset_index()

        if self.isForwardStrand:
            self.dfString["strand"] = self.dfString["min_pos"] % 3
        else:
            self.dfString["strand"] = self.dfString["max_pos"] % 3            

        self.addCodons()

        self.dfGeneString = pd.merge(self.dfString, self.geneTranscriptDf, on=['chromosome_identifier','gene_id', 'transcript_id'],how='left')
        self.dfGeneString["is_baseline"] = True if self.groupType == "baseline" else False
        self.dfGeneString["gene_predicted"] = False
        self.dfGeneString["same_strand"] = False
        self.dfGeneString = self.dfGeneString[self.reader.getGeneStringDfCols()]

    def writeGeneStringDf(self):
        self.dfGeneString = self.dfGeneString[self.reader.getGeneStringDfCols()]
        self.dfGeneString.to_csv(self.geneStringPath, index=False)

    def writeGeneStringCompleteDf(self):
        self.dfGeneString = self.dfGeneString[self.reader.getGeneStringDfCols()]
        self.dfGeneString.to_csv(self.geneStringCompletePath, index=False)

    def writeDf(self):
        self.df.to_csv(self.path, index=False)

    def getIntersectionGenes(self, commonGenes, getCommonValues):
        mask = self.dfGeneString["gene_string"].isin(commonGenes)
        return self.dfGeneString[mask if getCommonValues else ~mask].copy()

    def getGeneStringDf(self):
        return self.dfGeneString.copy()

    def updateMainDf(self, genePredictionDf):
        self.df = self.df.drop(columns=['predicted','gene_predicted'])
        self.df = pd.merge(self.df, genePredictionDf, on=['chromosome_identifier', 'gene_id', 'transcript_id', 'is_forward_strand'], how='left')
        self.writeDf()

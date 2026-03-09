import os
import pandas as pd
from sgevalviz.utils import *
import numpy as np

def defineFirstLastExon(df):
    isForward = df.iloc[0]["is_forward_strand"]
    g = df.groupby(["gene_id","transcript_id"])

    smaller = g.head(1).index
    bigger = g.tail(1).index

    firstExon = smaller if isForward else bigger
    lastExon = bigger if isForward else smaller

    df.loc[firstExon, "is_first_exon"] = True
    df.loc[lastExon, "is_last_exon"] = True

    return df

def defineIntronRetentionExon(df):
    df["exon_start_repeats"] = df.duplicated("region_start", keep=False)
    df["exon_end_repeats"] = df.duplicated("region_end", keep=False)

    min_end = df.groupby("region_start")["region_end"].transform("min")
    df["non_smallest_exon_from_same_start_group"] = (df["region_end"] != min_end)

    df["is_intron_retention_exon"] = (df["exon_start_repeats"] & df["exon_end_repeats"] & df["non_smallest_exon_from_same_start_group"])
    df.drop(columns=["exon_start_repeats", "exon_end_repeats", "non_smallest_exon_from_same_start_group"], inplace=True)

    return df

def dropLastIntron(df):
    dfLocal = df.groupby(['gene_id', 'transcript_id'])

    lastIntronId = dfLocal.tail(1).index
    df.drop(lastIntronId, inplace=True)

    return df

def unifiedIntronExonDf(dfExon,dfIntron):
    df = pd.concat([dfExon, dfIntron]).sort_index().reset_index(drop=True)
    df = df.sort_values(by=['gene_id', 'transcript_id', 'region_start']).reset_index(drop=True)

    df['next_exon_start'] = np.where(df['is_exon'], df['region_start'], np.nan)
    df['next_exon_start'] = df['next_exon_start'].bfill()
    df.loc[df['is_intron'], 'region_end'] = df.loc[df['is_intron'], 'next_exon_start'] - 1
    df.drop(columns='next_exon_start', inplace=True)

    return df

def getCodonDf(df,maskCol,regionStartValue):
    codonDf = df.loc[df[maskCol]]
    codonDf = codonDf[['chromosome_identifier','gene_id','transcript_id','region_start']]
    codonDf.rename(columns={'region_start':regionStartValue},inplace=True)

    return codonDf

def getDfString(df, startCodonDf, stopCodonDf):

    dfString = df.loc[df['is_exon']].copy()
    dfString['gene_string'] = dfString['region_start'].astype(str) + ';' + dfString['region_end'].astype(str)

    dfString = dfString.groupby(['chromosome_identifier','gene_id', 'transcript_id']
    ).agg(
        gene_string=('gene_string', '/'.join), # join all gene_string values
        exon_qtty=('gene_string', 'size'), # count how many were merged
        intron_retention_qtty=('is_intron_retention_exon','sum') #count how many of the exons are intron retention
    ).reset_index()

    dfString = pd.merge(dfString,startCodonDf,on=['chromosome_identifier','gene_id', 'transcript_id'],how='left')
    dfString = pd.merge(dfString,stopCodonDf,on=['chromosome_identifier','gene_id', 'transcript_id'],how='left')
    dfString[['start_codon_init', 'stop_codon_init', 'gene_string']] = (dfString[['start_codon_init', 'stop_codon_init', 'gene_string']].fillna(''))
    dfString['gene_string'] = "|" + dfString['start_codon_init'].astype(str) + '|' + dfString['gene_string'] + '|' + dfString['stop_codon_init'].astype(str) + "|"


    dfString['predicted'] = False
    dfString['gene_predicted'] = False

    return dfString

def getMasks(df,maskList):
    return [df[mask] for mask in maskList]


def getSubDfs(df,listOfPairs):
    dfsList = []
    for masks in listOfPairs:
        compositeMask = (masks[0] & masks[1]) if len(masks) == 2 else masks[0]
        newDf = df.loc[compositeMask].copy()
        dfsList.append(newDf)

    return dfsList

def fillCsv(saveFilesBasePath,csvPath,isBaseline,hasBothFiles,isForwardStrand):
    if not os.path.isfile(csvPath):
        return False

    df = pd.read_csv(csvPath)
    df = df.sort_values(by=['gene_id', 'transcript_id', 'region_start']).reset_index(drop=True)
    if isForwardStrand == None:
        isForwardStrand = df.iloc[0]["is_forward_strand"]

    # Initialize columns
    df['is_first_exon'] = False
    df['is_last_exon'] = False
    df['is_intron_retention_exon'] = False

    # Basic Masks
    isIntron, isExon = getMasks(df,["is_intron","is_exon"])
    isNotIntronOrExon = ~isIntron & ~isExon

    dfExon, dfIntron, dfNotIntronOrExon = getSubDfs(
        df,
        [
            [isExon],
            [isIntron],
            [isNotIntronOrExon]
        ]
    )

    dfExon = defineFirstLastExon(dfExon)
    dfExon = defineIntronRetentionExon(dfExon)
    dfIntron = dropLastIntron(dfIntron)

    dfExonIntron = unifiedIntronExonDf(dfExon, dfIntron)

    df = pd.concat([dfExonIntron, dfNotIntronOrExon]).sort_index().reset_index(drop=True)
    df['region_end'] = pd.to_numeric(df['region_end'], downcast='integer', errors='coerce')

    startCodonDf = getCodonDf(df,'is_start_codon','start_codon_init')
    stopCodonDf = getCodonDf(df,'is_stop_codon','stop_codon_init')

    dfString = getDfString(df,startCodonDf,stopCodonDf)

    df['predicted'] = False
    df['gene_predicted'] = False


    df.to_csv(csvPath, index=False)
    basePath = "/".join(csvPath.split('/')[:-1]) + "/"

    genePredictFilePath = f"{saveFilesBasePath}chromosomeCSVs/transcriptAndGeneBaselineFile.csv" if isBaseline else f"{saveFilesBasePath}chromosomeCSVs/transcriptAndGeneCandidateFile.csv" 
    genePredictDf = pd.read_csv(genePredictFilePath) 
    dfPredictString = pd.merge(dfString,genePredictDf,on=['chromosome_identifier','gene_id', 'transcript_id'],how='left')
    dfPredictString = dfPredictString[["chromosome_identifier","gene_id", "transcript_id","start_gene","end_gene","start_transcript","end_transcript","exon_qtty","intron_retention_qtty","gene_string","is_forward_strand","predicted"]]
    dfPredictString.to_csv(csvPath.split('.')[0]+"__gene_string.csv", index=False)
    if not hasBothFiles:
        dfPredictString["is_baseline"] = isBaseline
        dfPredictString["is_forward_strand"] = isForwardStrand
        dfPredictString["gene_predicted"] = False
        dfPredictString.to_csv(f"{basePath}/gene_transcript_predicted.csv", index=False)

    return df, dfPredictString, isForwardStrand

def fillCsvPredict(genePredictDf, csvDf, csvPath):
    genePredictDf = genePredictDf.drop(columns=['start_gene','end_gene','start_transcript','end_transcript','gene_string'])
    csvDf = csvDf.drop(columns=['predicted','gene_predicted'])

    fullCompariosonDf = pd.merge(csvDf,genePredictDf,on=['chromosome_identifier', 'gene_id', 'transcript_id', 'is_forward_strand'], how='left')

    fullCompariosonDf.to_csv(csvPath, index=False)

def filterDataframes(df,commonGenes,matchCommonValues,isBaseline):
    newDf = None
    mask = df['gene_string'].isin(commonGenes)
    newDf = df[mask if matchCommonValues else ~mask].copy()

    newDf["is_baseline"] = isBaseline

    return newDf

def predictedOrNotDf(candidateDf,baselineDf,predicted):
    newDf = pd.concat([candidateDf, baselineDf])
    newDf["predicted"] = predicted

    return newDf

def compareGenes(sf,baselineDf,candidateDf,isForwardStrand):

    # get unique values from both
    candidateGenes = candidateDf['gene_string'].unique()
    baselineGenes = baselineDf['gene_string'].unique()

    # intersection
    commonGenes = set(candidateGenes) & set(baselineGenes)

    # filter original DataFrames by baseline and common or not
    candidateCommon = filterDataframes(candidateDf,commonGenes,True,False)
    baselineCommon = filterDataframes(baselineDf,commonGenes,True,True)
    candidateNotCommon = filterDataframes(candidateDf,commonGenes,False,False)
    baselineNotCommon = filterDataframes(baselineDf,commonGenes,False,True)

    # predicted true/false
    predictedDf = predictedOrNotDf(candidateCommon, baselineCommon, True)
    notPredictedDf = predictedOrNotDf(candidateNotCommon, baselineNotCommon, False)


    # join predicted and non predicted
    fullDf = pd.concat([predictedDf,notPredictedDf]).copy()

    # get gene_predicted if any transcript predicts
    anyTranscriptPredictedDf = (
        fullDf.groupby('gene_id')['predicted']
        .any()
        .rename("gene_predicted")   # rename the Series itself
        .reset_index()              # turn it back into a DataFrame
    )
    fullDf = pd.merge(fullDf,anyTranscriptPredictedDf,on="gene_id", how="left")
    fullDf["is_forward_strand"] = isForwardStrand

    fullDf.to_csv(f"{sf}/gene_transcript_predicted.csv", index=False)

    return fullDf


def getChromosomeFolders(saveFilesBasePath):
    folder = f"{saveFilesBasePath}chromosomeCSVs"
    subfolders = [f.path for f in os.scandir(folder) if f.is_dir()]
    
    return subfolders

def hasBaseline(sf):
    return os.path.isfile(f"{sf}/processedBaselineFile.csv")

def hasCandidate(sf):
    return os.path.isfile(f"{sf}/processedCandidateFile.csv")

def appendFile(inputFile,outputFile):
    if not os.path.exists(inputFile):
        return

    with open(inputFile,"r") as f_in, open(outputFile,"a") as f_out:
        inputData = f_in.read().split("\n",1)[1]
        f_out.write(inputData)

def fillData(saveFilesBasePath,extraArgs):
    chromosomeFolders = getChromosomeFolders(saveFilesBasePath)
    for sf in chromosomeFolders:
        isForwardStrand = None
        if not (hasBaseline(sf) or hasCandidate(sf)):
            continue

        if hasBaseline(sf):
            baselineDf, baselineGeneStringDf, isForwardStrand = fillCsv(saveFilesBasePath,f"{sf}/processedBaselineFile.csv",True,hasCandidate(sf),isForwardStrand)

        if hasCandidate(sf):
            candidateDf, candidateGeneStringDf, isForwardStrand = fillCsv(saveFilesBasePath,f"{sf}/processedCandidateFile.csv",False,hasBaseline(sf),isForwardStrand)

        if hasBaseline(sf) and hasCandidate(sf):
            predictedDf = compareGenes(sf, baselineGeneStringDf, candidateGeneStringDf, isForwardStrand)
            fillCsvPredict(predictedDf,baselineDf,f"{sf}/processedBaselineFile.csv")
            fillCsvPredict(predictedDf,candidateDf,f"{sf}/processedCandidateFile.csv")


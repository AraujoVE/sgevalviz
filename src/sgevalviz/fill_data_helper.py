from __future__ import annotations
from sgevalviz.reader import Reader
from pandas import DataFrame, Series
import pandas as pd
import numpy as np
from itertools import product

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
        self.dfTranscript = None
        self.dfGenePredicted = None

    def getStrand(self):
        return "forward" if self.isForwardStrand else "reverse"

    def getDf(self):
        return self.df

    def getDfTranscript(self):
        return self.dfTranscript

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

        self.dfExon = self.dfExon.sort_values(by=['gene_id', 'transcript_id', 'region_start']).reset_index(drop=True)
        self.dfIntron = self.dfIntron.sort_values(by=['gene_id', 'transcript_id', 'region_start']).reset_index(drop=True)
        self.dfNotIntronOrExon = self.dfNotIntronOrExon.sort_values(by=['gene_id', 'transcript_id', 'region_start']).reset_index(drop=True)

    def defineFirstLastSingleExon(self):
        g = self.dfExon.groupby(["gene_id","transcript_id"])

        smaller = g.head(1).index
        bigger = g.tail(1).index

        firstExon = smaller if self.isForwardStrand else bigger
        lastExon = bigger if self.isForwardStrand else smaller

        self.dfExon.loc[firstExon, "is_first_exon"] = True
        self.dfExon.loc[lastExon, "is_last_exon"] = True
        self.dfExon["is_single_exon"] = self.dfExon["is_first_exon"] & self.dfExon["is_last_exon"]

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

    def setNucleotidesData(self):
        self.dfExon["nucleotide_size"] = (self.dfExon["region_end"] - self.dfExon["region_start"]) + 1
        self.dfIntron["nucleotide_size"] = (self.dfIntron["region_end"] - self.dfIntron["region_start"]) + 1

    def unifyDf(self):
        df = pd.concat([self.dfExon, self.dfIntron]).sort_index().reset_index(drop=True)
        df = df.sort_values(by=['gene_id', 'transcript_id', 'region_start']).reset_index(drop=True)

        df['next_exon_start'] = np.where(df['is_exon'], df['region_start'], np.nan)
        df['next_exon_start'] = df['next_exon_start'].bfill()
        df.loc[df['is_intron'], 'region_end'] = df.loc[df['is_intron'], 'next_exon_start'] - 1
        df.drop(columns='next_exon_start', inplace=True)

        self.df = pd.concat([df, self.dfNotIntronOrExon]).sort_index().reset_index(drop=True)
        self.df['region_end'] = pd.to_numeric(self.df['region_end'], downcast='integer', errors='coerce')

    def setTranscriptDf(self):
        self.df = self.df.sort_values(by=["gene_id", "transcript_id", "region_start"])


        df = self.df.assign(
            start_codon_pos=self.df["region_start"].where(self.df["is_start_codon"]),
            stop_codon_pos=self.df["region_start"].where(self.df["is_stop_codon"]),
            exon_start=self.df["region_start"].where(self.df["is_exon"]),
            exon_end=self.df["region_end"].where(self.df["is_exon"]),
            intron_start=self.df["region_start"].where(self.df["is_intron"]),
            intron_end=self.df["region_end"].where(self.df["is_intron"]),
            exon_size=self.df["nucleotide_size"].where(self.df["is_exon"]),
            intron_size=self.df["nucleotide_size"].where(self.df["is_intron"])
        )

        self.dfTranscript = (
            df.groupby(["gene_id", "transcript_id"])
            .agg(
                start_codon=("start_codon_pos", "min"),
                stop_codon=("stop_codon_pos", "min"),
                intron_starts=("intron_start", lambda s: s.dropna().tolist()),
                intron_ends=("intron_end", lambda s: s.dropna().tolist()),
                exon_starts=("exon_start", lambda s: s.dropna().tolist()),
                exon_ends=("exon_end", lambda s: s.dropna().tolist()),
                min_exon_start=("exon_start", "min"),
                max_exon_end=("exon_end", "max"),
                intron_retention_exons=("is_intron_retention_exon", "sum"),
                number_of_exons=("is_exon", "sum"),
                exon_avg_size=("exon_size", "mean"),
                intron_avg_size=("intron_size", "mean"),
                number_of_cds_nucleotides=("exon_size", "sum"),
                number_of_intron_nucleotides=("intron_size", "sum")
            )
            .reset_index()
        )
        #with pd.option_context("display.max_rows", None, "display.max_columns", None, "display.width", None, "display.expand_frame_repr", False):
        #    print("Main Df")
        #    print(self.df)
        #    print(self.dfTranscript)

        zipped = list(zip(self.dfTranscript["intron_starts"], self.dfTranscript["intron_ends"]))
        self.dfTranscript["introns"] = [
           list(zip(s, e)) for s, e in zip(self.dfTranscript["intron_starts"], self.dfTranscript["intron_ends"])
        ]

        self.dfTranscript["exons"] = [
           list(zip(s, e)) for s, e in zip(self.dfTranscript["exon_starts"], self.dfTranscript["exon_ends"])
        ]

        self.dfTranscript["cds_nucleotides"] = [
            set().union(*(range(int(s), int(e) + 1) for s, e in exons)) if exons else set()
            for exons in self.dfTranscript["exons"]
        ]

        exonList = self.dfTranscript["exons"]


        if self.isForwardStrand:
            self.dfTranscript["frame"] = self.dfTranscript["min_exon_start"] % 3
            firstExons = [l[0] if l else (None, None) for l in exonList]
            lastExons  = [l[-1] if l else (None, None) for l in exonList]
        else:
            self.dfTranscript["frame"] = self.dfTranscript["max_exon_end"] % 3
            firstExons = [l[-1] if l else (None, None) for l in exonList]
            lastExons  = [l[0] if l else (None, None) for l in exonList]

        self.dfTranscript["first_exon"] = firstExons
        self.dfTranscript["last_exon"] = lastExons

        self.dfTranscript["max_intron_retention"] = (
            self.dfTranscript.groupby("gene_id")["intron_retention_exons"]
            .transform("max")
        )

        self.dfTranscript["has_max_intron_retention"] = self.dfTranscript["max_intron_retention"] == self.dfTranscript["intron_retention_exons"]
        self.dfTranscript["gene_has_intron_retention"] = self.dfTranscript["max_intron_retention"] > 0

        self.dfTranscript["has_start_codon"] = self.dfTranscript["start_codon"].notna()
        self.dfTranscript["has_stop_codon"] = self.dfTranscript["stop_codon"].notna()


        self.dfTranscript.drop(["exon_starts", "exon_ends", "intron_retention_exons", "max_intron_retention"], axis=1, inplace=True)

        self.dfTranscript.rename(columns={
            "start_codon": f"{self.groupType}_start_codon",
            "has_start_codon": f"{self.groupType}_has_start_codon",
            "stop_codon": f"{self.groupType}_stop_codon",
            "has_stop_codon": f"{self.groupType}_has_stop_codon",
            "gene_id": f"{self.groupType}_gene_id",
            "transcript_id": f"{self.groupType}_transcript_id",
            "frame": f"{self.groupType}_frame",
            "introns": f"{self.groupType}_introns",
            "intron_starts": f"{self.groupType}_{'donnors' if self.isForwardStrand else 'acceptors'}",
            "intron_ends": f"{self.groupType}_{'acceptors' if self.isForwardStrand else 'donnors'}",
            "exons": f"{self.groupType}_exons",
            "first_exon": f"{self.groupType}_first_exon",
            "last_exon": f"{self.groupType}_last_exon",
            "has_max_intron_retention": f"{self.groupType}_has_max_intron_retention",
            "gene_has_intron_retention": f"{self.groupType}_gene_has_intron_retention",
            "number_of_exons": f"{self.groupType}_number_of_exons",
            "exon_avg_size": f"{self.groupType}_exon_avg_size",
            "intron_avg_size": f"{self.groupType}_intron_avg_size",
            "number_of_cds_nucleotides": f"{self.groupType}_number_of_cds_nucleotides",
            "number_of_intron_nucleotides": f"{self.groupType}_number_of_intron_nucleotides",
            "min_exon_start": f"{self.groupType}_cds_min",
            "max_exon_end": f"{self.groupType}_cds_max",
            "cds_nucleotides": f"{self.groupType}_cds_nucleotides"
        }, inplace=True)

    def getIntersectionSize(self, min1, max1, cds1, min2, max2, cds2):
        if not cds1 or not cds2:
            return 0
        if max1 < min2 or max2 < min1:
            return 0

        return len(cds1 & cds2)

    def getSetLen(self, listA, listB):
        return len(set(listA) & set(listB)) if listA and listB else 0

    def findCandidateTranscriptItsBaselineTranscript(self, baseline: FillDataHelper):
        baselineDf = baseline.getDfTranscript()
        candidateDf = self.dfTranscript

        candidates = list(zip(
            candidateDf["candidate_gene_id"],
            candidateDf["candidate_transcript_id"],
            candidateDf["candidate_cds_nucleotides"],
            candidateDf["candidate_cds_min"],
            candidateDf["candidate_cds_max"]
        ))
        
        baselines = list(zip(
            baselineDf["baseline_gene_id"],
            baselineDf["baseline_transcript_id"],
            baselineDf["baseline_cds_nucleotides"],
            baselineDf["baseline_cds_min"],
            baselineDf["baseline_cds_max"]
        ))

        results = {}

        for (cand_gene, cand_tx, cand_cds, cand_min, cand_max), (base_gene, base_tx, base_cds, base_min, base_max) in product(candidates, baselines):

            
            score = self.getIntersectionSize(cand_min, cand_max, cand_cds, base_min, base_max, base_cds)

            results_key = f"{cand_gene}___{cand_tx}"
            if (not results_key in results) or (results[results_key]["nucleotides_predicted"] < score):
                results[results_key] = {
                    "nucleotides_predicted": score,
                    "predicted": score > 0,
                    "totally_predicted": (score == len(cand_cds) and score == len(base_cds)),
                    "candidate_gene_id": cand_gene,
                    "baseline_gene_id": base_gene if score > 0 else None,
                    "candidate_transcript_id": cand_tx,
                    "baseline_transcript_id": base_tx if score > 0 else None,
                }

        baseDfGenePredicted = pd.DataFrame(list(results.values()))
        baseDfGenePredicted = baseDfGenePredicted.merge(candidateDf, on=["candidate_gene_id", "candidate_transcript_id"], how="left")
        self.dfGenePredicted = baseDfGenePredicted.merge(baselineDf, on=["baseline_gene_id", "baseline_transcript_id"], how="left")
        self.dfGenePredicted.drop(columns=["candidate_cds_nucleotides", "baseline_cds_nucleotides"], errors="ignore",  inplace=True)

        self.dfGenePredicted["start_codon_predicted"] = self.dfGenePredicted["candidate_start_codon"] == self.dfGenePredicted["baseline_start_codon"]
        self.dfGenePredicted["stop_codon_predicted"] = self.dfGenePredicted["candidate_stop_codon"] == self.dfGenePredicted["baseline_stop_codon"]
        self.dfGenePredicted["start_and_stop_codon_predicted"] = self.dfGenePredicted["start_codon_predicted"] & self.dfGenePredicted["stop_codon_predicted"]
        self.dfGenePredicted["first_exon_predicted"] = [
            fb == fc
            for fb, fc in zip(self.dfGenePredicted["baseline_first_exon"], self.dfGenePredicted["candidate_first_exon"]) 
        ]
        self.dfGenePredicted["last_exon_predicted"] = [
            lb == lc
            for lb, lc in zip(self.dfGenePredicted["baseline_last_exon"], self.dfGenePredicted["candidate_last_exon"]) 
        ]
        self.dfGenePredicted["first_and_last_exon_predicted"] = self.dfGenePredicted["first_exon_predicted"] & self.dfGenePredicted["last_exon_predicted"]

        for origin in ["baseline", "candidate"]:
            for subtype in ["donnors", "acceptors", "exons", "introns"]:
                col = f"{origin}_{subtype}"
                self.dfGenePredicted[col] = [
                    x if isinstance(x, list) else [] for x in self.dfGenePredicted[col]
                ]



        self.dfGenePredicted["donnors_predicted"] = [
            self.getSetLen(db, dc)
            for db, dc in zip(self.dfGenePredicted["baseline_donnors"], self.dfGenePredicted["candidate_donnors"])
        ]
        self.dfGenePredicted["acceptors_predicted"] = [
            self.getSetLen(ab, ac)
            for ab, ac in zip(self.dfGenePredicted["baseline_acceptors"], self.dfGenePredicted["candidate_acceptors"])
        ]
        self.dfGenePredicted["exons_predicted"] = [
            self.getSetLen(eb, ec)
            for eb, ec in zip(self.dfGenePredicted["baseline_exons"], self.dfGenePredicted["candidate_exons"])
        ]
        self.dfGenePredicted["introns_predicted"] = [
            self.getSetLen(ib, ic)
            for ib, ic in zip(self.dfGenePredicted["baseline_introns"], self.dfGenePredicted["candidate_introns"])
        ]

        return self.dfGenePredicted


    def enrinchDf(self):
        self.setIntronExonAndOthers()
        self.defineFirstLastSingleExon()
        self.defineIntronRetentionExon()
        self.dropLastIntron()
        self.setNucleotidesData()
        self.unifyDf()
        self.setTranscriptDf()

    def writeDf(self):
        self.df.to_csv(self.path, index=False)

    def getGeneStringDf(self):
        return self.dfGeneString.copy()

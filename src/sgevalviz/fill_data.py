import pandas as pd
from sgevalviz.reader import Reader
from sgevalviz.fill_data_helper import FillDataHelper

def boolTo01List(boolDf):
    return [1 if i == True else 0 for i in boolDf.tolist()]

def fillCsv(reader: Reader, chromosomePath, groupType, hasBothFiles):

    dfHelper = FillDataHelper(chromosomePath, groupType, reader)

    dfHelper.enrinchDf()
    dfHelper.writeDf()

    return dfHelper

def divideByTotallyAndPartiallyPredicted(strand, reader: Reader, predictedDf, baselineDfHelper: FillDataHelper, region):
    baselineTranscriptsDf = baselineDfHelper.getDfTranscript().copy()
    baselineExonQttyDf = baselineTranscriptsDf.groupby("baseline_gene_id", as_index=False
    ).agg(
        baseline_min_exon_qtty=("baseline_number_of_exons","min"),
        baseline_max_exon_qtty=("baseline_number_of_exons","max"),
        baseline_transcripts_qtty=("baseline_number_of_exons","size")
    )

    copyDf = pd.merge(predictedDf, baselineExonQttyDf, on=["baseline_gene_id"], how="left")
    totallyPredictedDf = copyDf[copyDf["totally_predicted"]].copy()
    partiallyPredictedDf = copyDf[~copyDf["totally_predicted"]].copy()

    dfs = {
        "partially_predicted_gene": partiallyPredictedDf,
        "totally_predicted_gene": totallyPredictedDf
    }

    sameStartFrame = (partiallyPredictedDf["candidate_start_frame"] == partiallyPredictedDf["baseline_start_frame"])
    sameStartFrameList = boolTo01List(sameStartFrame)

    sameEndFrame = (partiallyPredictedDf["candidate_end_frame"] == partiallyPredictedDf["baseline_end_frame"])
    sameEndFrameList = boolTo01List(sameEndFrame)

    sameStartAndEndFrameList = boolTo01List(sameStartFrame & sameEndFrame)

    exonRatiosList = (partiallyPredictedDf["candidate_number_of_exons"] / partiallyPredictedDf["baseline_number_of_exons"]).tolist()
    nucleotideRatiosList = (partiallyPredictedDf["candidate_number_of_cds_nucleotides"] / partiallyPredictedDf["baseline_number_of_cds_nucleotides"]).tolist()

    for curStrand in [strand, "general"]:
        reader.statisticsUpdate([curStrand, "correct_frame_on_first_nucleotide"], sameStartFrameList, True, region)
        reader.statisticsUpdate([curStrand, "correct_frame_on_last_nucleotide"], sameEndFrameList, True, region)
        reader.statisticsUpdate([curStrand, "correct_frame_on_first_and_last_nucleotide"], sameStartAndEndFrameList, True, region)

        reader.statisticsUpdate([curStrand, "ratio_of_number_of_exons_in_prediction_per_reference"], exonRatiosList, False, region)
        reader.statisticsUpdate([curStrand, "ratio_of_number_of_nucleotides_in_prediction_per_reference"], nucleotideRatiosList, False, region)

    for predType, df in dfs.items():
        intronRetentionList = boolTo01List(df["baseline_gene_has_intron_retention"])
        maxIntronRetentionList = boolTo01List(df["baseline_gene_has_intron_retention"] & df["baseline_has_max_intron_retention"])

        hasLeastExonsInVariants = ((df["baseline_min_exon_qtty"] == df["baseline_number_of_exons"]) & (df["baseline_transcripts_qtty"] > 1))
        hasLeastExonsList = boolTo01List(hasLeastExonsInVariants)

        hasMostExonsInVariants = ((df["baseline_max_exon_qtty"] == df["baseline_number_of_exons"]) & (df["baseline_transcripts_qtty"] > 1))
        hasMostExonsList = boolTo01List(hasMostExonsInVariants)

        numberOfExonsList = (df["baseline_number_of_exons"]).tolist()
        sizeOfExonsList = (df["baseline_exon_avg_size"]).tolist()
        sizeOfIntronsList = (df["baseline_intron_avg_size"]).tolist()

        for curStrand in [strand, "general"]:
            for pred in [predType, "predicted_gene"]:
                reader.statisticsUpdate([curStrand, pred, "prediction_has_intron_retention"], intronRetentionList, True, region)
                reader.statisticsUpdate([curStrand, pred, "prediction_maximizes_intron_retention"], maxIntronRetentionList, True, region)

                reader.statisticsUpdate([curStrand, pred, "splice_variant_selected_has_the_least_number_of_cds"], hasLeastExonsList, True, region)
                reader.statisticsUpdate([curStrand, pred, "splice_variant_selected_has_the_most_number_of_cds"], hasMostExonsList, True, region)

                reader.statisticsUpdate([curStrand, pred, "number_of_cds_in_splice_variant_selected"], numberOfExonsList, False, region)
                reader.statisticsUpdate([curStrand, pred, "average_cds_size_in_splice_variant_selected"], sizeOfExonsList, False, region)
                reader.statisticsUpdate([curStrand, pred, "average_intron_size_in_splice_variant_selected"], sizeOfIntronsList, False, region)


def getUnpredictedDf(baseDf, unpredictedGeneId, isCandidate):
    if len(unpredictedGeneId) == 0:
        return None

    baseName = "candidate" if isCandidate else "baseline"

    unpredictedDf = baseDf[baseDf[f"{baseName}_gene_id"].isin(unpredictedGeneId)].copy()
    unpredictedDf = unpredictedDf.sort_values(by=f"{baseName}_number_of_cds_nucleotides").drop_duplicates(subset=[f"{baseName}_gene_id"], keep="first")

    unpredictedDf["nucleotides_predicted"] = 0
    unpredictedDf["predicted"] = False
    unpredictedDf["totally_predicted"] = False
    unpredictedDf["start_codon_predicted"] = False 
    unpredictedDf["stop_codon_predicted"] = False 
    unpredictedDf["start_and_stop_codon_predicted"] = False 
    unpredictedDf["first_exon_predicted"] = False 
    unpredictedDf["last_exon_predicted"] = False 
    unpredictedDf["first_and_last_exon_predicted"] = False
    unpredictedDf["donnors_predicted"] = 0 
    unpredictedDf["acceptors_predicted"] = 0 
    unpredictedDf["exons_predicted"] = 0 
    unpredictedDf["introns_predicted"] = 0

    unpredictedDf = unpredictedDf[[
        f"{baseName}_gene_id", f"{baseName}_transcript_id", "totally_predicted",
        "nucleotides_predicted", f"{baseName}_number_of_cds_nucleotides",
        "start_codon_predicted", f"{baseName}_has_start_codon", "stop_codon_predicted", f"{baseName}_has_stop_codon", "start_and_stop_codon_predicted",
        "introns_predicted", "exons_predicted", f"{baseName}_number_of_exons",
        "first_exon_predicted", "last_exon_predicted", "first_and_last_exon_predicted",
        "donnors_predicted", "acceptors_predicted"
    ]].copy()

    return unpredictedDf

def getPredictedDf(dfPrediction, unpredictedDf, isCandidate):
    baseName = "candidate" if isCandidate else "baseline"
    predictedData = dfPrediction[[
        f"{baseName}_gene_id", f"{baseName}_transcript_id", "totally_predicted",
        "nucleotides_predicted", f"{baseName}_number_of_cds_nucleotides",
        "start_codon_predicted", f"{baseName}_has_start_codon", "stop_codon_predicted", f"{baseName}_has_stop_codon", "start_and_stop_codon_predicted",
        "introns_predicted", "exons_predicted", f"{baseName}_number_of_exons",
        "first_exon_predicted", "last_exon_predicted", "first_and_last_exon_predicted",
        "donnors_predicted", "acceptors_predicted"
    ]].copy()

    concatenatedDf = predictedData if unpredictedDf is None else pd.concat([predictedData, unpredictedDf], ignore_index=True)

    concatenatedDf = concatenatedDf.sort_values(
        by="nucleotides_predicted",
        ascending=False
    ).drop_duplicates(
        subset=[f"{baseName}_gene_id"], 
        keep="first"
    )

    return concatenatedDf

def noMultiExonStatistics(df, reader: Reader, starterText, strand, recOrPre, singleOrMultipleString, region):
    totallyPredictedGenes = df.groupby(f"{starterText}_gene_id", as_index=False)["totally_predicted"].any()

    totallyPredictedGenesList = boolTo01List(totallyPredictedGenes["totally_predicted"])
    totalNucleotidesList = boolTo01List(df[f"{starterText}_number_of_cds_nucleotides"])
    totalStartCodonList = boolTo01List(df[f"{starterText}_has_start_codon"])
    totalStopCodonList = boolTo01List(df[f"{starterText}_has_stop_codon"])
    totalStartAndStopCodonList = boolTo01List(df[f"{starterText}_has_start_codon"] & df[f"{starterText}_has_stop_codon"])

    for curStrand in [strand, "general"]:
        for singleOrMultiple in [singleOrMultipleString, "single_or_multiple_exon"]:
            reader.statisticsUpdate([curStrand, singleOrMultiple, recOrPre, "gene"], totallyPredictedGenesList, True, region)
            reader.statisticsUpdate([curStrand, singleOrMultiple, recOrPre, "nucleotide"], totalNucleotidesList, True, region)
            reader.statisticsUpdate([curStrand, singleOrMultiple, recOrPre, "start_codon"], totalStartCodonList, True, region)
            reader.statisticsUpdate([curStrand, singleOrMultiple, recOrPre, "stop_codon"], totalStopCodonList, True, region)
            reader.statisticsUpdate([curStrand, singleOrMultiple, recOrPre, "start_and_stop_codon"], totalStartAndStopCodonList, True, region)


def noMultiExonStatistics(df, reader: Reader, starterText, strand, recOrPre, singleOrMultipleString, region):
    #Gene Recall/Precision takes into account only per gene_id
    totallyPredictedGenes = df.groupby(f"{starterText}_gene_id", as_index=False)["totally_predicted"].any()
    totallyPredictedGenesList = boolTo01List(totallyPredictedGenes["totally_predicted"])

    #Nucleotides/Exons Recall/Precision only takes into account predicted ones
    hasPredictionDf = df[df["predicted"]]
    nucleotidesBool = (hasPredictionDf["nucleotides_predicted"] / hasPredictionDf[f"{starterText}_number_of_cds_nucleotides"]).tolist()
    exonsBool = (hasPredictionDf["exons_predicted"] / hasPredictionDf[f"{starterText}_number_of_exons"]).tolist()


    #Start/Stop Codon Recall/Precision only takes into account when there IS the start/stop codon prediction
    startCodonBool = boolTo01List(df.loc[df[f"{starterText}_has_start_codon"], "start_codon_predicted"])
    stopCodonBool = boolTo01List(df.loc[df[f"{starterText}_has_stop_codon"], "stop_codon_predicted"]) 

    startAndStopCodonHasStart = df.loc[(df[f"{starterText}_has_start_codon"]) & (df[f"{starterText}_has_stop_codon"]), "start_codon_predicted"]
    startAndStopCodonHasStop = df.loc[(df[f"{starterText}_has_start_codon"]) & (df[f"{starterText}_has_stop_codon"]), "stop_codon_predicted"]
    startAndStopCodonBool = boolTo01List(startAndStopCodonHasStart & startAndStopCodonHasStop)

    for curStrand in [strand, "general"]:
        for singleOrMultiple in [singleOrMultipleString, "single_or_multiple_exon"]:
            reader.statisticsUpdate([curStrand, singleOrMultiple, recOrPre, "gene"], totallyPredictedGenesList, True, region)
            reader.statisticsUpdate([curStrand, singleOrMultiple, recOrPre, "nucleotide"], nucleotidesBool, True, region)
            reader.statisticsUpdate([curStrand, singleOrMultiple, recOrPre, "exon"], exonsBool, True, region)
            reader.statisticsUpdate([curStrand, singleOrMultiple, recOrPre, "start_codon"], startCodonBool, True, region)
            reader.statisticsUpdate([curStrand, singleOrMultiple, recOrPre, "stop_codon"], stopCodonBool, True, region)
            reader.statisticsUpdate([curStrand, singleOrMultiple, recOrPre, "start_and_stop_codon"], startAndStopCodonBool, True, region)




def onlyMultiExonStatistics(df, reader: Reader, starterText, strand, recOrPre, region):
    hasPredictionDf = df[df["predicted"]].copy()

    #Recall/Precision only takes into account predicted ones
    intronsBool = (hasPredictionDf["introns_predicted"] / hasPredictionDf[f"{starterText}_number_of_introns"]).tolist()
    firstExonsBool = boolTo01List(hasPredictionDf["first_exon_predicted"])
    lastExonsBool = boolTo01List(hasPredictionDf["last_exon_predicted"])
    firstAndLastExonsBool = boolTo01List(hasPredictionDf["first_exon_predicted"] & hasPredictionDf["last_exon_predicted"])
    donnorsBool = (hasPredictionDf["donnors_predicted"] / hasPredictionDf[f"{starterText}_number_of_introns"]).tolist()
    acceptorsBool = (hasPredictionDf["acceptors_predicted"] / hasPredictionDf[f"{starterText}_number_of_introns"]).tolist()

    for curStrand in [strand, "general"]:
        for singleOrMultiple in ["multiple_exon", "single_or_multiple_exon"]:
            reader.statisticsUpdate([curStrand, singleOrMultiple, recOrPre, "intron"], intronsBool, True, region)
            reader.statisticsUpdate([curStrand, singleOrMultiple, recOrPre, "first_exon"], firstExonsBool, True, region)
            reader.statisticsUpdate([curStrand, singleOrMultiple, recOrPre, "last_exon"], lastExonsBool, True, region)
            reader.statisticsUpdate([curStrand, singleOrMultiple, recOrPre, "first_and_last_exon"], firstAndLastExonsBool, True, region)
            reader.statisticsUpdate([curStrand, singleOrMultiple, recOrPre, "donnors"], donnorsBool, True, region)
            reader.statisticsUpdate([curStrand, singleOrMultiple, recOrPre, "acceptors"], acceptorsBool, True, region)

def recallOrPrecision(df, reader: Reader, strand, recOrPre, region):
    candidateOrBaseline = "baseline" if recOrPre == "gene_predicted_recall" else "candidate"
    singleExonOnlyDf = df[df[f"{candidateOrBaseline}_number_of_exons"] == 1]
    multipleExonOnlyDf = df[df[f"{candidateOrBaseline}_number_of_exons"] > 1]

    totalLen = len(df)
    singleExonLen = len(singleExonOnlyDf)
    multiExonLen = totalLen - singleExonLen
    singleExonList = [0]*multiExonLen + [1]*singleExonLen

    for curStrand in [strand, "general"]:
        reader.statisticsUpdate([curStrand, recOrPre, "single_exon_occurance"], singleExonList, True, region)

    if len(singleExonOnlyDf) > 0:
        noMultiExonStatistics(singleExonOnlyDf, reader, candidateOrBaseline, strand, recOrPre, "single_exon", region)
    if len(multipleExonOnlyDf) > 0:
        noMultiExonStatistics(multipleExonOnlyDf, reader, candidateOrBaseline, strand, recOrPre, "multiple_exon", region)
        onlyMultiExonStatistics(multipleExonOnlyDf, reader, candidateOrBaseline, strand, recOrPre, region)

    return

def recallAndPrecision(strand, reader: Reader, candidateDfHelper: FillDataHelper, baselineDfHelper: FillDataHelper, dfPrediction, region):
    uniqueCandidate = candidateDfHelper.getDfTranscript()
    uniqueBaseline = baselineDfHelper.getDfTranscript()

    allCandidateGeneIds = set(uniqueCandidate["candidate_gene_id"])
    predictedCandidateGeneIds = set(dfPrediction["candidate_gene_id"])
    unpredictedCandidateGeneIds = allCandidateGeneIds - predictedCandidateGeneIds

    allBaselineGeneIds = set(uniqueBaseline["baseline_gene_id"])
    predictedBaselineGeneIds = set(dfPrediction["baseline_gene_id"])
    unpredictedBaselineGeneIds = allBaselineGeneIds - predictedBaselineGeneIds


    unpredictedCandidates = getUnpredictedDf(uniqueCandidate, unpredictedCandidateGeneIds, True)
    precisionDf = getPredictedDf(dfPrediction, unpredictedCandidates, True)

    unpredictedBaselines = getUnpredictedDf(uniqueBaseline, unpredictedBaselineGeneIds, False)
    recallDf = getPredictedDf(dfPrediction, unpredictedBaselines, False)

    recallOrPrecision(recallDf, reader, strand, "recall_of", region)
    recallOrPrecision(precisionDf, reader, strand, "precision_of", region)

    return

def findPrediction(reader: Reader, candidateDfHelper: FillDataHelper, baselineDfHelper: FillDataHelper, dfPrediction, region):
    strand = candidateDfHelper.getStrand()

    #statisticsUpdate__Strand__Value
    allBaselineGenes = len(set(baselineDfHelper.getDf()["gene_id"].dropna().unique()))
    predictedGenes = len(set(dfPrediction["baseline_gene_id"].dropna().unique()))
    unpredictedBaselineGenes = allBaselineGenes - predictedGenes
    predictionList = [0]*unpredictedBaselineGenes + [1]*predictedGenes

    for curStrand in [strand, "general"]:
        reader.statisticsUpdate([curStrand, "genes_ignored"], predictionList, True, region)
        reader.statisticsUpdate([curStrand, "prediction_on_wrong_strand"], [0], True, region)

    dfPrediction = dfPrediction[dfPrediction["predicted"]]
    #statisticsUpdate__StrandRefGenePartPred__Value
    divideByTotallyAndPartiallyPredicted(strand, reader, dfPrediction, baselineDfHelper, region)
    recallAndPrecision(strand, reader, candidateDfHelper, baselineDfHelper, dfPrediction, region)

    return

def addEmptyData(reader: Reader, baselineDfHelper: FillDataHelper, candidateDfHelper: FillDataHelper, hasBaseline, region):
    strand = baselineDfHelper.getStrand() if hasBaseline else candidateDfHelper.getStrand()

    if hasBaseline:
        allBaselineGenes = len(set(baselineDfHelper.getDf()["gene_id"].dropna().unique()))
        for curStrand in [strand, "general"]:
            reader.statisticsUpdate([curStrand, "genes_ignored"], [1]*allBaselineGenes, True, region)
    else:
        for curStrand in [strand, "general"]:
            reader.statisticsUpdate([curStrand, "prediction_on_wrong_strand"], [1], True, region)


def fillData(reader: Reader):
    chromosomeFolders = reader.getChromosomeFoldersList()
    for region, chromosomePath in enumerate(chromosomeFolders):
        hasCandidate, hasBaseline = reader.hasGtfFile(chromosomePath, "candidate"), reader.hasGtfFile(chromosomePath, "baseline")
        hasBothFiles = hasCandidate and hasBaseline
        if not (hasBaseline or hasCandidate):
            continue

        if hasBaseline:
            baselineDfHelper = fillCsv(reader, chromosomePath, "baseline", hasBothFiles)

        if hasCandidate:
            candidateDfHelper = fillCsv(reader, chromosomePath, "candidate", hasBothFiles)

        if hasBothFiles:
            dfPrediction = candidateDfHelper.findCandidateTranscriptItsBaselineTranscript(baselineDfHelper)
            findPrediction(reader, candidateDfHelper, baselineDfHelper, dfPrediction, region)
        else:
            addEmptyData(reader, baselineDfHelper, candidateDfHelper, hasBaseline, region)

    #reader.setFinalResults()

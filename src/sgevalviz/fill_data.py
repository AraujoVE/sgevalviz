import pandas as pd
from sgevalviz.reader import Reader
from sgevalviz.fill_data_helper import FillDataHelper

def boolTo01List(boolDf):
    return [1 if i == True else 0 for i in boolDf.tolist()]

def valueDividendAndDivisor(boolDf):
    divisor = len(boolDf)
    dividend = boolDf.sum()
    value = 0 if divisor == 0 else dividend / divisor

    return value, dividend, divisor

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
    sameStartFrameValue, sameStartFrameDividend, sameStartFrameDivisor = valueDividendAndDivisor(sameStartFrame)

    sameEndFrame = (partiallyPredictedDf["candidate_end_frame"] == partiallyPredictedDf["baseline_end_frame"])
    sameEndFrameValue, sameEndFrameDividend, sameEndFrameDivisor = valueDividendAndDivisor(sameEndFrame)

    sameStartAndEndFrameValue, sameStartAndEndFrameDividend, sameStartAndEndFrameDivisor = valueDividendAndDivisor(sameStartFrame & sameEndFrame)

    exonCandidateList, exonBaselineList, exonRatiosList = partiallyPredictedDf["candidate_number_of_exons"].tolist(), partiallyPredictedDf["baseline_number_of_exons"].tolist(), (partiallyPredictedDf["candidate_number_of_exons"] / partiallyPredictedDf["baseline_number_of_exons"]).tolist()
    nucleotideCandidateList, nucleotideBaselineList, nucleotideRatiosList = partiallyPredictedDf["candidate_number_of_cds_nucleotides"].tolist(), partiallyPredictedDf["baseline_number_of_cds_nucleotides"].tolist(), (partiallyPredictedDf["candidate_number_of_cds_nucleotides"] / partiallyPredictedDf["baseline_number_of_cds_nucleotides"]).tolist()

    for curStrand in [strand, "general"]:
        reader.statisticsUpdate([curStrand, "correct_frame_on_first_nucleotide"], region, True, True, False, sameStartFrameValue, sameStartFrameDividend, sameStartFrameDivisor)
        reader.statisticsUpdate([curStrand, "correct_frame_on_last_nucleotide"], region, True, True, False, sameEndFrameValue, sameEndFrameDividend, sameEndFrameDivisor)
        reader.statisticsUpdate([curStrand, "correct_frame_on_first_and_last_nucleotide"], region, True, True, False, sameStartAndEndFrameValue, sameStartAndEndFrameDividend, sameStartAndEndFrameDivisor)

        reader.statisticsUpdate([curStrand, "ratio_of_number_of_exons_in_prediction_per_reference"], region, False, True, True, exonRatiosList, exonCandidateList, exonBaselineList)
        reader.statisticsUpdate([curStrand, "ratio_of_number_of_nucleotides_in_prediction_per_reference"], region, False, True, True, nucleotideRatiosList, nucleotideCandidateList, nucleotideBaselineList)

    for predType, df in dfs.items():
        hasIntronRetentionDf = df[df["baseline_gene_has_intron_retention"]].copy()
        if len(hasIntronRetentionDf) == 0:
            intronRetentionValue, intronRetentionDividend, intronRetentionDivisor = 0, 0, 0
            maxIntronRetentionValue, maxIntronRetentionDividend, maxIntronRetentionDivisor = 0, 0, 0
        else:
            intronRetentionValue, intronRetentionDividend, intronRetentionDivisor = valueDividendAndDivisor(hasIntronRetentionDf["baseline_has_intron_retention"])
            maxIntronRetentionValue, maxIntronRetentionDividend, maxIntronRetentionDivisor = valueDividendAndDivisor(hasIntronRetentionDf["baseline_has_max_intron_retention"])

        hasLeastExonsInVariants = ((df["baseline_min_exon_qtty"] == df["baseline_number_of_exons"]) & (df["baseline_transcripts_qtty"] > 1))
        hasLeastExonsValue, hasLeastExonsDividend, hasLeastExonsDivisor = valueDividendAndDivisor(hasLeastExonsInVariants)

        hasMostExonsInVariants = ((df["baseline_max_exon_qtty"] == df["baseline_number_of_exons"]) & (df["baseline_transcripts_qtty"] > 1))
        hasMostExonsValue, hasMostExonsDividend, hasMostExonsDivisor = valueDividendAndDivisor(hasMostExonsInVariants)

        numberOfExonsList = (df["baseline_number_of_exons"]).tolist()
        sizeOfExonsList = (df["baseline_exon_avg_size"]).tolist()
        sizeOfIntronsList = (df["baseline_intron_avg_size"].dropna()).tolist()

        for curStrand in [strand, "general"]:
            for pred in [predType, "predicted_gene"]:
                reader.statisticsUpdate([curStrand, pred, "prediction_has_intron_retention"], region, True, True, False, intronRetentionValue, intronRetentionDividend, intronRetentionDivisor)
                reader.statisticsUpdate([curStrand, pred, "prediction_maximizes_intron_retention"], region, True, True, False, maxIntronRetentionValue, maxIntronRetentionDividend, maxIntronRetentionDivisor)

                reader.statisticsUpdate([curStrand, pred, "splice_variant_selected_has_the_least_number_of_cds"], region, True, True, False, hasLeastExonsValue, hasLeastExonsDividend, hasLeastExonsDivisor)
                reader.statisticsUpdate([curStrand, pred, "splice_variant_selected_has_the_most_number_of_cds"], region, True, True, False, hasMostExonsValue, hasMostExonsDividend, hasMostExonsDivisor)

                reader.statisticsUpdate([curStrand, pred, "number_of_cds_in_splice_variant_selected"], region, False, False, True, numberOfExonsList, [], [])
                reader.statisticsUpdate([curStrand, pred, "average_cds_size_in_splice_variant_selected"], region, False, False, True, sizeOfExonsList, [], [])
                reader.statisticsUpdate([curStrand, pred, "average_intron_size_in_splice_variant_selected"], region, False, False, True, sizeOfIntronsList, [], [])


def getUnpredictedDf(baseDf, unpredictedGeneId, isCandidate):
    if len(unpredictedGeneId) == 0:
        return None

    baseName = "candidate" if isCandidate else "baseline"

    unpredictedDf = baseDf[baseDf[f"{baseName}_gene_id"].isin(unpredictedGeneId)].copy()
    unpredictedDf = unpredictedDf.sort_values(by=f"{baseName}_number_of_cds_nucleotides").drop_duplicates(subset=[f"{baseName}_gene_id"], keep="first")

    unpredictedDf["nucleotides_predicted"] = 0
    unpredictedDf["predicted"] = False
    unpredictedDf["totally_predicted"] = False
    unpredictedDf["predicted"] = False
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
        f"{baseName}_gene_id", f"{baseName}_transcript_id", "totally_predicted", "predicted",
        "nucleotides_predicted", f"{baseName}_number_of_cds_nucleotides",
        "start_codon_predicted", f"{baseName}_has_start_codon", "stop_codon_predicted", f"{baseName}_has_stop_codon", "start_and_stop_codon_predicted",
        "introns_predicted", "exons_predicted", f"{baseName}_number_of_exons", f"{baseName}_number_of_introns",
        "first_exon_predicted", "last_exon_predicted", "first_and_last_exon_predicted",
        "donnors_predicted", "acceptors_predicted"
    ]].copy()

    return unpredictedDf

def getPredictedDf(dfPrediction, unpredictedDf, isCandidate):
    baseName = "candidate" if isCandidate else "baseline"
    predictedData = dfPrediction[[
        f"{baseName}_gene_id", f"{baseName}_transcript_id", "totally_predicted", "predicted",
        "nucleotides_predicted", f"{baseName}_number_of_cds_nucleotides",
        "start_codon_predicted", f"{baseName}_has_start_codon", "stop_codon_predicted", f"{baseName}_has_stop_codon", "start_and_stop_codon_predicted",
        "introns_predicted", "exons_predicted", f"{baseName}_number_of_exons",  f"{baseName}_number_of_introns",
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
    #Gene Recall/Precision takes into account only per gene_id
    totallyPredictedGenes = df.groupby(f"{starterText}_gene_id", as_index=False)["totally_predicted"].any()
    totallyPredictedGenesValue, totallyPredictedGenesDividend, totallyPredictedGenesDivisor = valueDividendAndDivisor(totallyPredictedGenes["totally_predicted"])

    #Nucleotides/Exons Recall/Precision only takes into account predicted ones
    hasPredictionDf = df[df["predicted"]]

    nucleotidesBool = (hasPredictionDf["nucleotides_predicted"] / hasPredictionDf[f"{starterText}_number_of_cds_nucleotides"]).tolist()
    nucleotidesDividend = hasPredictionDf["nucleotides_predicted"].tolist()
    nucleotidesDivisor = hasPredictionDf[f"{starterText}_number_of_cds_nucleotides"].tolist()

    exonsBool = (hasPredictionDf["exons_predicted"] / hasPredictionDf[f"{starterText}_number_of_exons"]).tolist()
    exonsDividend = hasPredictionDf["exons_predicted"].tolist()
    exonsDivisor = hasPredictionDf[f"{starterText}_number_of_exons"].tolist()


    #Start/Stop Codon Recall/Precision only takes into account when there IS the start/stop codon prediction
    startCodonBoolValue, startCodonBoolDividend, startCodonBoolDivisor = valueDividendAndDivisor(df.loc[df[f"{starterText}_has_start_codon"], "start_codon_predicted"])
    stopCodonBoolValue, stopCodonBoolDividend, stopCodonBoolDivisor = valueDividendAndDivisor(df.loc[df[f"{starterText}_has_stop_codon"], "stop_codon_predicted"]) 

    startAndStopCodonHasStart = df.loc[(df[f"{starterText}_has_start_codon"]) & (df[f"{starterText}_has_stop_codon"]), "start_codon_predicted"]
    startAndStopCodonHasStop = df.loc[(df[f"{starterText}_has_start_codon"]) & (df[f"{starterText}_has_stop_codon"]), "stop_codon_predicted"]
    startAndStopCodonBoolValue, startAndStopCodonBoolDividend, startAndStopCodonBoolDivisor = valueDividendAndDivisor(startAndStopCodonHasStart & startAndStopCodonHasStop)

    for curStrand in [strand, "general"]:
        for singleOrMultiple in [singleOrMultipleString, "single_or_multiple_exon"]:
            reader.statisticsUpdate([curStrand, singleOrMultiple, recOrPre, "gene"], region, True, True, False, totallyPredictedGenesValue, totallyPredictedGenesDividend, totallyPredictedGenesDivisor)
            reader.statisticsUpdate([curStrand, singleOrMultiple, recOrPre, "nucleotide"], region, True, True, True, nucleotidesBool, nucleotidesDividend, nucleotidesDivisor)
            reader.statisticsUpdate([curStrand, singleOrMultiple, recOrPre, "exon"], region, True, True, True, exonsBool, exonsDividend, exonsDivisor)
            reader.statisticsUpdate([curStrand, singleOrMultiple, recOrPre, "start_codon"], region, True, True, False, startCodonBoolValue, startCodonBoolDividend, startCodonBoolDivisor)
            reader.statisticsUpdate([curStrand, singleOrMultiple, recOrPre, "stop_codon"], region, True, True, False, stopCodonBoolValue, stopCodonBoolDividend, stopCodonBoolDivisor)
            reader.statisticsUpdate([curStrand, singleOrMultiple, recOrPre, "start_and_stop_codon"], region, True, True, False, startAndStopCodonBoolValue, startAndStopCodonBoolDividend, startAndStopCodonBoolDivisor)

def onlyMultiExonStatistics(df, reader: Reader, starterText, strand, recOrPre, region):
    hasPredictionDf = df[df["predicted"]].copy()

    #Recall/Precision only takes into account predicted ones
    intronsBool = (hasPredictionDf["introns_predicted"] / hasPredictionDf[f"{starterText}_number_of_introns"]).tolist()
    intronsDividend = hasPredictionDf["introns_predicted"].tolist()

    firstExonsBoolValue, firstExonsBoolDividend, firstExonsBoolDivisor = valueDividendAndDivisor(hasPredictionDf["first_exon_predicted"])
    lastExonsBoolValue, lastExonsBoolDividend, lastExonsBoolDivisor = valueDividendAndDivisor(hasPredictionDf["last_exon_predicted"])
    firstAndLastExonsBoolValue, firstAndLastExonsBoolDividend, firstAndLastExonsBoolDivisor = valueDividendAndDivisor(hasPredictionDf["first_exon_predicted"] & hasPredictionDf["last_exon_predicted"])

    donnorsBool = (hasPredictionDf["donnors_predicted"] / hasPredictionDf[f"{starterText}_number_of_introns"]).tolist()
    donnorsDividend = hasPredictionDf["donnors_predicted"].tolist()

    acceptorsBool = (hasPredictionDf["acceptors_predicted"] / hasPredictionDf[f"{starterText}_number_of_introns"]).tolist()
    acceptorsDividend = hasPredictionDf["acceptors_predicted"].tolist()

    intronsDonnorAcceptorDivisor = hasPredictionDf[f"{starterText}_number_of_introns"].tolist()

    for curStrand in [strand, "general"]:
        for singleOrMultiple in ["multiple_exon", "single_or_multiple_exon"]:
            reader.statisticsUpdate([curStrand, singleOrMultiple, recOrPre, "intron"], region, True, True, True, intronsBool, intronsDividend, intronsDonnorAcceptorDivisor)
            reader.statisticsUpdate([curStrand, singleOrMultiple, recOrPre, "first_exon"], region, True, True, False, firstExonsBoolValue, firstExonsBoolDividend, firstExonsBoolDivisor)
            reader.statisticsUpdate([curStrand, singleOrMultiple, recOrPre, "last_exon"], region, True, True, False, lastExonsBoolValue, lastExonsBoolDividend, lastExonsBoolDivisor)
            reader.statisticsUpdate([curStrand, singleOrMultiple, recOrPre, "first_and_last_exon"], region, True, True, False, firstAndLastExonsBoolValue, firstAndLastExonsBoolDividend, firstAndLastExonsBoolDivisor)
            reader.statisticsUpdate([curStrand, singleOrMultiple, recOrPre, "donnors"], region, True, True, True, donnorsBool, donnorsDividend, intronsDonnorAcceptorDivisor)
            reader.statisticsUpdate([curStrand, singleOrMultiple, recOrPre, "acceptors"], region, True, True, True, acceptorsBool, acceptorsDividend, intronsDonnorAcceptorDivisor)

def recallOrPrecision(df, reader: Reader, strand, recOrPre, region):
    candidateOrBaseline = "baseline" if recOrPre == "recall_of" else "candidate"
    singleExonOnlyDf = df[df[f"{candidateOrBaseline}_number_of_exons"] == 1]
    multipleExonOnlyDf = df[df[f"{candidateOrBaseline}_number_of_exons"] > 1]

    totalLen = len(df)
    singleExonLen = len(singleExonOnlyDf)
    singleExonPerc = singleExonLen / totalLen

    for curStrand in [strand, "general"]:
        reader.statisticsUpdate([curStrand, recOrPre, "single_exon_occurance"], region, True, True, False, singleExonPerc, singleExonLen, totalLen)

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
    prediction = unpredictedBaselineGenes / allBaselineGenes

    for curStrand in [strand, "general"]:
        reader.statisticsUpdate([curStrand, "genes_ignored"], region, True, True, False, prediction, unpredictedBaselineGenes, allBaselineGenes)
        reader.statisticsUpdate([curStrand, "no_reference_on_prediction_region"], region, True, True, False, 0, 0, 1)
        reader.statisticsUpdate([curStrand, "no_prediction_on_reference_region"], region, True, True, False, 0, 0, 1)

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
            reader.statisticsUpdate([curStrand, "genes_ignored"], region, True, True, False, allBaselineGenes, allBaselineGenes, allBaselineGenes)
            reader.statisticsUpdate([curStrand, "no_prediction_on_reference_region"], region, True, True, False, 1, 1, 1)
    else:
        for curStrand in [strand, "general"]:
            reader.statisticsUpdate([curStrand, "no_reference_on_prediction_region"], region, True, True, False, 1, 1, 1)


def fillData(reader: Reader):
    chromosomeFolders = reader.getChromosomeFoldersList()
    for region, chromosomePath in enumerate(chromosomeFolders):
        hasCandidate, hasBaseline = reader.hasGtfFile(chromosomePath, "candidate"), reader.hasGtfFile(chromosomePath, "baseline")
        hasBothFiles = hasCandidate and hasBaseline
        baselineDfHelper, candidateDfHelper = None, None
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

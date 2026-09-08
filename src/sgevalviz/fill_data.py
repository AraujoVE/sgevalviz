import pandas as pd
from sgevalviz.reader import Reader
from sgevalviz.fill_data_helper import FillDataHelper

def fillCsv(reader: Reader, chromosomePath, groupType, hasBothFiles):

    dfHelper = FillDataHelper(chromosomePath, groupType, reader)

    dfHelper.enrinchDf()
    dfHelper.writeDf()

    return dfHelper

def divideByTotallyAndPartiallyPredicted(strand, reader: Reader, predictedDf, baselineDfHelper: FillDataHelper):
    totallyPredictedDf = predictedDf[predictedDf["totally_predicted"]].copy()
    partiallyPredictedDf = predictedDf[~predictedDf["totally_predicted"]].copy()
    partiallyPredictedDfSize = len(partiallyPredictedDf)

    dfs = {
        "reference_gene_partially_predicted": partiallyPredictedDf,
        "reference_gene_totally_predicted": totallyPredictedDf
    }


    sameFrameCount = (partiallyPredictedDf["candidate_frame"] == partiallyPredictedDf["baseline_frame"]).sum()
    exonRatiosSummed = (partiallyPredictedDf["candidate_number_of_exons"] / partiallyPredictedDf["baseline_number_of_exons"]).sum()
    nucleotideRatiosSummed = (partiallyPredictedDf["candidate_number_of_cds_nucleotides"] / partiallyPredictedDf["baseline_number_of_cds_nucleotides"]).sum()

    for curStrand in [strand, "general"]:
        reader.statisticsUpdate__StrandRefGenePartPred__Value(curStrand, "selected_model_transcript_on_same_frame_as_selected_reference_transcript-percentage", sameFrameCount, partiallyPredictedDfSize)
        reader.statisticsUpdate__StrandRefGenePartPred__Value(curStrand, "ratio_of_number_of_exons_in_model_selected_transcript_per_reference_selected_transcript-average", exonRatiosSummed, partiallyPredictedDfSize)
        reader.statisticsUpdate__StrandRefGenePartPred__Value(curStrand, "ratio_of_number_of_nucleotides_in_model_selected_transcript_per_reference_selected_transcript-average", nucleotideRatiosSummed, partiallyPredictedDfSize)

    baselineTranscriptsDf = baselineDfHelper.getDfTranscript()

    for predType, df in dfs.items():
        hasIntronRetention = df["baseline_gene_has_intron_retention"].sum()
        hasMaxIntronRetention = df["baseline_has_max_intron_retention"].sum()

        for curStrand in [strand, "general"]:
            for pred in [predType, "reference_gene_predicted"]:
                reader.statisticsUpdate__StrandRefGenePredHasIntronRetExonInModel__Value(curStrand, pred, hasMaxIntronRetention, hasIntronRetention)

        selectedGenes = df["baseline_gene_id"].unique()
        selectedTranscripts = df[["baseline_gene_id", "baseline_transcript_id"]].copy()
        selectedTranscripts["transcript_predicted"] = True

        curBaselineDf = baselineTranscriptsDf[baselineTranscriptsDf["baseline_gene_id"].isin(selectedGenes)].copy()
        curBaselineDf = pd.merge(curBaselineDf, selectedTranscripts, on=["baseline_gene_id", "baseline_transcript_id"], how="left")
        curBaselineDf["transcript_predicted"] = curBaselineDf["transcript_predicted"].fillna(False)

        dfSelected = curBaselineDf[curBaselineDf["transcript_predicted"] == True].groupby("baseline_gene_id"
        ).agg(
            avg_exon_qtty=("baseline_number_of_exons","mean"),
            avg_exon_size=("baseline_exon_avg_size","mean"),
            avg_intron_size=("baseline_intron_avg_size","mean")
        )

        dfGeneral = curBaselineDf.groupby("baseline_gene_id"
        ).agg(
            avg_exon_qtty=("baseline_number_of_exons","mean"),
            avg_exon_size=("baseline_exon_avg_size","mean"),
            avg_intron_size=("baseline_intron_avg_size","mean")
        )

        for data in [{"df": dfSelected, "string": "selected_model_transcript"}, {"df": dfGeneral, "string": "average_model_transcript"}]:
            curDf = data["df"]
            modelStr = data["string"]
            exonQtty = curDf["avg_exon_qtty"].sum()
            exonSize = curDf["avg_exon_size"].sum()
            intronSize = curDf["avg_intron_size"].sum()
            totLen = len(curDf)

            for curStrand in [strand, "general"]:
                for curPred in [predType, "reference_gene_predicted"]:
                    reader.statisticsUpdate__Strand_RefGenePred_ModelTransc__Value(curStrand, curPred, modelStr, "average_size_of_exons_per_transcript-average", exonSize, totLen)
                    reader.statisticsUpdate__Strand_RefGenePred_ModelTransc__Value(curStrand, curPred, modelStr, "average_size_of_introns_per_transcript-average", intronSize, totLen)
                    reader.statisticsUpdate__Strand_RefGenePred_ModelTransc__Value(curStrand, curPred, modelStr, "number_of_exons_per_transcript-average", exonQtty, totLen)


def getUnpredictedDf(predictedPairs, baseDf, candidateOrBaselinePairs, isCandidate):
    baseName = "candidate" if isCandidate else "baseline"

    unpredictedDf = baseDf[[pair not in predictedPairs for pair in candidateOrBaselinePairs]].copy()
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

    concatenatedDf = pd.concat([predictedData, unpredictedDf], ignore_index=True)

    concatenatedDf = concatenatedDf.sort_values(
        by="nucleotides_predicted",
        ascending=False
    ).drop_duplicates(
        subset=[f"{baseName}_gene_id", f"{baseName}_transcript_id"], 
        keep="first"
    )

    return concatenatedDf

def noMultiExonStatistics(df, reader: Reader, starterText, strand, recOrPre, singleOrMultipleString):
    totalGenes = len(df[f"{starterText}_gene_id"].unique())
    totallyPredictedGenes = df.groupby(f"{starterText}_gene_id")["totally_predicted"].any().sum()

    totalNucleotides = df[f"{starterText}_number_of_cds_nucleotides"].sum()
    predictedNucleotides = df["nucleotides_predicted"].sum()

    totalStartCodon = df[f"{starterText}_has_start_codon"].sum()
    predictedStartCodon = df["start_codon_predicted"].sum()

    totalStopCodon = df[f"{starterText}_has_stop_codon"].sum()
    predictedStopCodon = df["stop_codon_predicted"].sum()

    totalStartAndStopCodon = (df[f"{starterText}_has_start_codon"] & df[f"{starterText}_has_stop_codon"]).sum() 
    predictedStartAndStopCodon = df["start_and_stop_codon_predicted"].sum()

    for curStrand in [strand, "general"]:
        for singleOrMultiple in [singleOrMultipleString, "any_quantity_exon_selected_model_transcript"]:
            reader.statisticsUpdate__Strand_RefGenePredRecallPrecision_ExonNumberPerTranscript__Value(curStrand, recOrPre, singleOrMultiple,"totally_predicted_genes_prediction-percentage", totallyPredictedGenes, totalGenes)
            reader.statisticsUpdate__Strand_RefGenePredRecallPrecision_ExonNumberPerTranscript__Value(curStrand, recOrPre, singleOrMultiple, "nucleotide_prediction-percentage", predictedNucleotides, totalNucleotides)
            reader.statisticsUpdate__Strand_RefGenePredRecallPrecision_ExonNumberPerTranscript__Value(curStrand, recOrPre, singleOrMultiple, "start_codon_prediction-percentage", predictedStartCodon, totalStartCodon)
            reader.statisticsUpdate__Strand_RefGenePredRecallPrecision_ExonNumberPerTranscript__Value(curStrand, recOrPre, singleOrMultiple, "stop_codon_prediction-percentage", predictedStopCodon, totalStopCodon)
            reader.statisticsUpdate__Strand_RefGenePredRecallPrecision_ExonNumberPerTranscript__Value(curStrand, recOrPre, singleOrMultiple, "same_time_start_codon_and_stop_codon_prediction-percentage", predictedStartAndStopCodon, totalStartAndStopCodon)



def onlyMultiExonStatistics(df, reader: Reader, starterText, strand, recOrPre):
    numberOfRows = len(df)
    numberOfExons = df[f"{starterText}_number_of_exons"].sum()
    numberOfIntrons = numberOfExons - numberOfRows

    intronsPredicted = df["introns_predicted"].sum()
    exonsPredicted = df["exons_predicted"].sum()
    firstExonsPredicted = df["first_exon_predicted"].sum()
    lastExonsPredicted = df["last_exon_predicted"].sum()
    firstAndLastExonsPredicted = df["first_and_last_exon_predicted"].sum()
    donnorsPredicted = df["donnors_predicted"].sum()
    acceptorsPredicted = df["acceptors_predicted"].sum()

    for curStrand in [strand, "general"]:
        for multipleOrAny in ["multiple_exon_selected_model_transcript", "any_quantity_exon_selected_model_transcript"]:
            reader.statisticsUpdate__Strand_RefGenePredRecallPrecision_ExonNumberPerTranscript__Value(curStrand, recOrPre, multipleOrAny,"average_intron_prediction-percentage", intronsPredicted, numberOfIntrons)
            reader.statisticsUpdate__Strand_RefGenePredRecallPrecision_ExonNumberPerTranscript__Value(curStrand, recOrPre, multipleOrAny,"average_exon_prediction-percentage", exonsPredicted, numberOfExons)
            reader.statisticsUpdate__Strand_RefGenePredRecallPrecision_ExonNumberPerTranscript__Value(curStrand, recOrPre, multipleOrAny,"first_exon_prediction-percentage", firstExonsPredicted, numberOfRows)
            reader.statisticsUpdate__Strand_RefGenePredRecallPrecision_ExonNumberPerTranscript__Value(curStrand, recOrPre, multipleOrAny,"last_exon_prediction-percentage", lastExonsPredicted, numberOfRows)
            reader.statisticsUpdate__Strand_RefGenePredRecallPrecision_ExonNumberPerTranscript__Value(curStrand, recOrPre, multipleOrAny,"same_time_first_exon_and_last_exon_prediction-percentage", firstAndLastExonsPredicted, numberOfRows)
            reader.statisticsUpdate__Strand_RefGenePredRecallPrecision_ExonNumberPerTranscript__Value(curStrand, recOrPre, multipleOrAny,"average_donnor_prediction-percentage", donnorsPredicted, numberOfIntrons)
            reader.statisticsUpdate__Strand_RefGenePredRecallPrecision_ExonNumberPerTranscript__Value(curStrand, recOrPre, multipleOrAny,"average_acceptor_prediction-percentage", acceptorsPredicted, numberOfIntrons)



def recallOrPrecision(df, reader: Reader, strand, recOrPre):
    candidateOrBaseline = "baseline" if recOrPre == "gene_predicted_recall" else "candidate"
    singleExonOnlyDf = df[df[f"{candidateOrBaseline}_number_of_exons"] == 1]
    multipleExonOnlyDf = df[df[f"{candidateOrBaseline}_number_of_exons"] > 1]

    totalLen = len(df)
    singleExonLen = len(singleExonOnlyDf)

    for curStrand in [strand, "general"]:
        reader.statisticsUpdate__StrandRefGenePredRecallPrecision__Value(curStrand, recOrPre, singleExonLen, totalLen)

    if len(singleExonOnlyDf) > 0:
        noMultiExonStatistics(singleExonOnlyDf, reader, candidateOrBaseline, strand, recOrPre, "single_exon_selected_model_transcript")
    if len(multipleExonOnlyDf) > 0:
        noMultiExonStatistics(multipleExonOnlyDf, reader, candidateOrBaseline, strand, recOrPre, "multiple_exon_selected_model_transcript")
        onlyMultiExonStatistics(multipleExonOnlyDf, reader, candidateOrBaseline, strand, recOrPre)

    return

def recallAndPrecision(strand, reader: Reader, candidateDfHelper: FillDataHelper, baselineDfHelper: FillDataHelper, dfPrediction):
    predictedPairsCandidate = set(zip(dfPrediction["candidate_gene_id"], dfPrediction["candidate_transcript_id"]))
    predictedPairsBaseline = set(zip(dfPrediction["baseline_gene_id"], dfPrediction["baseline_transcript_id"]))

    uniqueCandidate = candidateDfHelper.getDfTranscript()
    candidatePairs = list(zip(uniqueCandidate["candidate_gene_id"], uniqueCandidate["candidate_transcript_id"]))

    uniqueBaseline = baselineDfHelper.getDfTranscript()
    baselinePairs = list(zip(uniqueBaseline["baseline_gene_id"], uniqueBaseline["baseline_transcript_id"]))

    unpredictedCandidates = getUnpredictedDf(predictedPairsCandidate, uniqueCandidate, candidatePairs, True)
    precisionDf = getPredictedDf(dfPrediction, unpredictedCandidates, True)

    unpredictedBaselines = getUnpredictedDf(predictedPairsBaseline, uniqueBaseline, baselinePairs, False)
    recallDf = getPredictedDf(dfPrediction, unpredictedBaselines, False)

    recallOrPrecision(recallDf, reader, strand, "gene_predicted_recall")
    recallOrPrecision(precisionDf, reader, strand, "gene_predicted_precision")

    return

def findPrediction(reader: Reader, candidateDfHelper: FillDataHelper, baselineDfHelper: FillDataHelper, dfPrediction):
    strand = candidateDfHelper.getStrand()

    #statisticsUpdate__Strand__Value
    allBaselineGenes = len(set(baselineDfHelper.getDf()["gene_id"].dropna().unique()))
    predictedGenes = len(set(dfPrediction["baseline_gene_id"].dropna().unique()))
    unpredictedBaselineGenes = allBaselineGenes - predictedGenes

    for curStrand in [strand, "general"]:
        reader.statisticsUpdate__Strand__Value(curStrand, "reference_gene_unpredicted-percentage", unpredictedBaselineGenes, allBaselineGenes)
        reader.statisticsUpdate__Strand__Value(curStrand, "no_prediction_for_reference_gene_on_same_strand-percentage", 0, 1)

    dfPrediction = dfPrediction[dfPrediction["predicted"]]
    #statisticsUpdate__StrandRefGenePartPred__Value
    divideByTotallyAndPartiallyPredicted(strand, reader, dfPrediction, baselineDfHelper)
    recallAndPrecision(strand, reader, candidateDfHelper, baselineDfHelper, dfPrediction)

    return

def addEmptyData(reader: Reader, baselineDfHelper: FillDataHelper):
    strand = baselineDfHelper.getStrand()

    allBaselineGenes = len(set(baselineDfHelper.getDf()["gene_id"].dropna().unique()))
    for curStrand in [strand, "general"]:
        reader.statisticsUpdate__Strand__Value(curStrand, "reference_gene_unpredicted-percentage", allBaselineGenes, allBaselineGenes)
        reader.statisticsUpdate__Strand__Value(curStrand, "no_prediction_for_reference_gene_on_same_strand-percentage", 1, 1)
    

def fillData(reader: Reader):
    chromosomeFolders = reader.getChromosomeFoldersList()
    for chromosomePath in chromosomeFolders:
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
            findPrediction(reader, candidateDfHelper, baselineDfHelper, dfPrediction)
        elif hasBaseline:
            addEmptyData(reader, baselineDfHelper)

    reader.setFinalResults()
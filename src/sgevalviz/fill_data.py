import pandas as pd
from sgevalviz.reader import Reader
from sgevalviz.fill_data_helper import FillDataHelper

def fillCsv(reader: Reader, chromosomePath, groupType, hasBothFiles):

    dfHelper = FillDataHelper(chromosomePath, groupType, reader)

    dfHelper.enrinchDf()
    dfHelper.writeDf()

    dfHelper.generateGeneStringDf()
    dfHelper.writeGeneStringDf()

    #TODO: Probably wrong path
    if not hasBothFiles:
        dfHelper.writeGeneStringCompleteDf()

    return dfHelper

def predictedOrNotDf(candidateDf,baselineDf,predicted):
    newDf = pd.concat([candidateDf, baselineDf])
    newDf["predicted"] = predicted
    newDf.drop(columns="gene_predicted", inplace=True)

    return newDf

def checkIfValuesCross(a_start, a_end, b_start, b_end):
    return not (a_end < b_start or b_end < a_start)

def addSameStrandDf(candidateDfHelper: FillDataHelper, baselineDfHelper: FillDataHelper, genePredictionDf):
    candidateDf, baselineDf = (
        candidateDfHelper.getGeneStringDf(),
        baselineDfHelper.getGeneStringDf(),
    )

    candidateDf["strand"] = candidateDf["strand"].astype(int)
    baselineDf["strand"] = baselineDf["strand"].astype(int)


    candidateStrands = [
        candidateDf.loc[candidateDf["strand"] == i, ["min_pos", "max_pos"]].to_numpy()
        for i in range(3)
    ]

    baselineDfUnpredicted = baselineDf.loc[~baselineDf["gene_predicted"]].copy()

    baselineDfUnpredicted["same_strand"] = baselineDfUnpredicted.apply(
        lambda row: any(
            checkIfValuesCross(row["min_pos"], row["max_pos"], minP, maxP)
            for minP, maxP in candidateStrands[row["strand"]]
        ),
        axis=1,
    )

    sameStrandDf = (
        baselineDfUnpredicted
        .groupby("gene_id", as_index=False)["same_strand"]
        .any()
    )

    baselineDfPredicted = (
        baselineDf.loc[baselineDf["gene_predicted"], ["gene_id"]]
        .drop_duplicates()
        .assign(same_strand=True)
    )

    candidateGenes = (
        candidateDf[["gene_id"]]
        .drop_duplicates()
        .assign(same_strand=False)
    )

    sameStrandDf = (
        pd.concat(
            [baselineDfPredicted, sameStrandDf, candidateGenes],
            ignore_index=True,
        )
        .drop_duplicates("gene_id", keep="first")
    )

    genePredictionDf.drop(columns="same_strand", inplace=True)
    genePredictionDf = pd.merge(genePredictionDf, sameStrandDf, on='gene_id', how='left')

    return genePredictionDf

def getGenePrediction(candidateDfHelper: FillDataHelper, baselineDfHelper: FillDataHelper, commonGenes):
    candidateCommon = candidateDfHelper.getIntersectionGenes(commonGenes, True)
    baselineCommon = baselineDfHelper.getIntersectionGenes(commonGenes, True)
    candidateNotCommon = candidateDfHelper.getIntersectionGenes(commonGenes, False)
    baselineNotCommon = baselineDfHelper.getIntersectionGenes(commonGenes, False)

    predictedDf = predictedOrNotDf(candidateCommon, baselineCommon, True)
    notPredictedDf = predictedOrNotDf(candidateNotCommon, baselineNotCommon, False)

    genePredictionDf = pd.concat([predictedDf,notPredictedDf]).copy()
    anyTranscriptPredictedDf = (
        genePredictionDf.groupby('gene_id')['predicted']
        .any()
        .rename("gene_predicted")   # rename the Series itself
        .reset_index()              # turn it back into a DataFrame
    )
    genePredictionDf = pd.merge(genePredictionDf,anyTranscriptPredictedDf,on="gene_id", how="left")

    return genePredictionDf


def findPrediction(reader: Reader, candidateDfHelper: FillDataHelper, baselineDfHelper: FillDataHelper, chromosomePath):
    candidateGenes, baselineGenes = candidateDfHelper.getUniqueGeneString(), baselineDfHelper.getUniqueGeneString()
    commonGenes = set(candidateGenes) & set(baselineGenes)
    genePredictionDf = getGenePrediction(candidateDfHelper, baselineDfHelper, commonGenes)
    genePredictionDf = addSameStrandDf(candidateDfHelper, baselineDfHelper, genePredictionDf)
    genePredictionDf = genePredictionDf[reader.getGeneStringDfCols()]

    genePredictionDf.to_csv(reader.getDefinedChromosomeSingleGeneStringPath(chromosomePath), index=False)

    genePredictionDf = genePredictionDf[['chromosome_identifier', 'gene_id', 'transcript_id', 'is_forward_strand', 'predicted', 'gene_predicted']]
    
    candidateDfHelper.updateMainDf(genePredictionDf)
    baselineDfHelper.updateMainDf(genePredictionDf)

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
            findPrediction(reader, candidateDfHelper, baselineDfHelper, chromosomePath)

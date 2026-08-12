import pandas as pd
import json
import os
from sgevalviz.utils import *
import itertools
import numpy as np
from sgevalviz.reader import Reader

def getChromosomeFolders(saveFilesBasePath):
    folder = f"{saveFilesBasePath}chromosomeCSVs"
    subfolders = [f.path for f in os.scandir(folder) if f.is_dir()]
    
    return subfolders

def getMasks(df,maskStrList):
    return [df[col].to_numpy() for col in maskStrList]

def groupDf(df,groupbyList,aggTerm):
    return df.groupby(groupbyList)[aggTerm].any().reset_index()

def addSizeParams(dfList):
    for df in dfList:
        df["region_length"] = (df["region_end"] - df["region_start"]).abs() + 1

def intersectMasks(mainMask, secondaryMaskList):
    return [(mainMask & secondaryMask) for secondaryMask in secondaryMaskList]

def maskPairs(listOfMaskPairs):
    return [(m1 & m2) for m1, m2 in listOfMaskPairs]

def countTrueValues(masks):
    return [mask.sum() for mask in masks]

def countTotalSizeSum(regionLength, maskList):
    return [(regionLength[mask]).sum() for mask in maskList]

def getStrandedVersionOfList(preDfListGeneral,isForwardStrand):
    strandString = "forward" if isForwardStrand == True else "reverse"
    preDfListStranded = [(f"{i}_{strandString}",v) for i, v in preDfListGeneral]

    return preDfListStranded

def getEmptyDf(hasForwardStrand, hasReverseStrand, hasTotal, useDf):
    basePredictionData = [(pair,0) for pair in [
        "exons_total",
        "introns_total",
        "start_codons_total",
        "stop_codons_total",
        "first_exons_only_total",
        "last_exons_only_total",
        "single_exons_total",
        "donors",
        "acceptors",
        "exons_predicted",
        "exons_unpredicted",
        "introns_predicted",
        "introns_unpredicted",
        "start_codons_predicted",
        "start_codons_unpredicted",
        "stop_codons_predicted",
        "stop_codons_unpredicted",
        "first_exons_only_predicted",
        "first_exons_only_unpredicted",
        "last_exons_only_predicted",
        "last_exons_only_unpredicted",
        "single_exons_predicted",
        "single_exons_unpredicted",
        "donors_predicted",
        "donors_unpredicted",
        "acceptors_predicted",
        "acceptors_unpredicted",
        "exons_total_size_sum",
        "introns_total_size_sum",
        "first_exons_only_total_size_sum",
        "last_exons_only_total_size_sum",
        "single_exons_total_size_sum",
        "exons_predicted_size_sum",
        "exons_unpredicted_size_sum",
        "introns_predicted_size_sum",
        "introns_unpredicted_size_sum",
        "first_exons_only_predicted_size_sum",
        "first_exons_only_unpredicted_size_sum",
        "last_exons_only_predicted_size_sum",
        "last_exons_only_unpredicted_size_sum",
        "single_exons_predicted_size_sum",
        "single_exons_unpredicted_size_sum"
    ]]
    predictionData = basePredictionData if hasTotal else []
    predictionData += getStrandedVersionOfList(basePredictionData,True) if hasForwardStrand else []
    predictionData += getStrandedVersionOfList(basePredictionData,False) if hasReverseStrand else []

    predictionDataDf = pd.DataFrame(predictionData, columns=["identifier", "value"])

    return predictionDataDf if useDf else predictionData 

def regionPredictionData(mainDf, isForwardStrand):
    predictionData = []

    # Intermediate Dfs
    exonDf = mainDf.loc[mainDf['is_exon']]
    intronDf = mainDf.loc[mainDf['is_intron']]

    # Grouped Dfs
    mainDfDedup = groupDf(mainDf,["is_exon","is_intron","is_start_codon","is_stop_codon","region_start","region_end"],"region_predicted")
    firstsDf = groupDf(exonDf,["is_first_exon","is_last_exon","region_start","region_end"],"region_predicted")
    donorDf = groupDf(intronDf,["region_start"],"donor_predicted")
    acceptorDf = groupDf(intronDf,["region_end"],"acceptor_predicted")

    #Add size params
    addSizeParams([mainDfDedup,firstsDf])

    # Region Length
    mainDfDedupRegionLength = mainDfDedup["region_length"].to_numpy()
    firstsDfRegionLength = firstsDf["region_length"].to_numpy()

    # Base Masks
    firstDfPredicted, isFirstExon, isLastExon = getMasks(firstsDf,["region_predicted","is_first_exon","is_last_exon"]) 
    mainRegionPredicted, isExon, isIntron, isStartCodon, isStopCodon = getMasks(mainDfDedup,["region_predicted","is_exon","is_intron","is_start_codon","is_stop_codon"])
    isDonorPredicted = getMasks(donorDf,["donor_predicted"])[0]
    isAcceptorPredicted = getMasks(acceptorDf,["acceptor_predicted"])[0]

    # Secondary Masks
    isNotFirstExon, isNotLastExon = ~isFirstExon, ~isLastExon
    isFirstExonOnly, isLastExonOnly, isSingleExon = maskPairs([[isFirstExon,isNotLastExon],[isNotFirstExon,isLastExon],[isFirstExon,isLastExon]])
    isExonPredicted, isIntronPredicted, isStartCodonPredicted, isStopCodonPredicted = intersectMasks(mainRegionPredicted,[isExon,isIntron,isStartCodon,isStopCodon])
    isFirstExonOnlyPredicted, isLastExonOnlyPredicted, isSingleExonPredicted = intersectMasks(firstDfPredicted,[isFirstExonOnly,isLastExonOnly,isSingleExon])

    # Count Quantities
    totalDonors, totalAcceptors = len(donorDf), len(acceptorDf)
    totalExons, totalIntrons, totalStartCodons, totalStopCodons, totalFirstExonsOnly, totalLastExonsOnly, totalSingleExons, successExons, successIntrons, successStartCodons, successStopCodons, successFirstExonsOnly, successLastExonsOnly, successSingleExons, successDonors, successAcceptors = countTrueValues([isExon,isIntron,isStartCodon,isStopCodon,isFirstExonOnly,isLastExonOnly,isSingleExon,isExonPredicted,isIntronPredicted,isStartCodonPredicted,isStopCodonPredicted,isFirstExonOnlyPredicted,isLastExonOnlyPredicted,isSingleExonPredicted,isDonorPredicted,isAcceptorPredicted])

    # Count Sizes
    totalExonsSizeSum, totalIntronsSizeSum, successExonsSizeSum, successIntronsSizeSum = countTotalSizeSum(mainDfDedupRegionLength,[isExon,isIntron,isExonPredicted,isIntronPredicted])
    totalFirstExonsOnlySizeSum, totalLastExonsOnlySizeSum, totalSingleExonsSizeSum, successFirstExonsOnlySizeSum, successLastExonsOnlySizeSum, successSingleExonsSizeSum = countTotalSizeSum(firstsDfRegionLength,[isFirstExonOnly,isLastExonOnly,isSingleExon,isFirstExonOnlyPredicted,isLastExonOnlyPredicted,isSingleExonPredicted])

    # Write Data
    predictionDataGeneral = [
        #Total Count of Values, in general
        ("exons_total", totalExons),
        ("introns_total", totalIntrons),
        ("start_codons_total", totalStartCodons),
        ("stop_codons_total", totalStopCodons),
        ("first_exons_only_total", totalFirstExonsOnly),
        ("last_exons_only_total", totalLastExonsOnly),
        ("single_exons_total", totalSingleExons),
        ("donors_total", totalDonors),
        ("acceptors_total", totalAcceptors),        
        #Total Count of Values, Predicted and Unpredicted
        ("exons_predicted", successExons),
        ("exons_unpredicted", totalExons - successExons),
        ("introns_predicted", successIntrons),
        ("introns_unpredicted", totalIntrons - successIntrons),
        ("start_codons_predicted", successStartCodons),
        ("start_codons_unpredicted", totalStartCodons - successStartCodons),
        ("stop_codons_predicted", successStopCodons),
        ("stop_codons_unpredicted", totalStopCodons - successStopCodons),
        ("first_exons_only_predicted", successFirstExonsOnly),
        ("first_exons_only_unpredicted", totalFirstExonsOnly - successFirstExonsOnly),
        ("last_exons_only_predicted", successLastExonsOnly),
        ("last_exons_only_unpredicted", totalLastExonsOnly - successLastExonsOnly),
        ("single_exons_predicted", successSingleExons),
        ("single_exons_unpredicted", totalSingleExons - successSingleExons),
        ("donors_predicted", successDonors),
        ("donors_unpredicted", totalDonors - successDonors),
        ("acceptors_predicted", successAcceptors),
        ("acceptors_unpredicted", totalAcceptors - successAcceptors),
        #Size sum of Values, in general
        ("exons_total_size_sum", totalExonsSizeSum),
        ("introns_total_size_sum", totalIntronsSizeSum),
        ("first_exons_only_total_size_sum", totalFirstExonsOnlySizeSum),
        ("last_exons_only_total_size_sum", totalLastExonsOnlySizeSum),
        ("single_exons_total_size_sum", totalSingleExonsSizeSum),
        #Size sum of Values, Predicted and Unpredicted
        ("exons_predicted_size_sum", successExonsSizeSum),
        ("exons_unpredicted_size_sum", totalExonsSizeSum - successExonsSizeSum),
        ("introns_predicted_size_sum", successIntronsSizeSum),
        ("introns_unpredicted_size_sum", totalIntronsSizeSum - successIntronsSizeSum),
        ("first_exons_only_predicted_size_sum", successFirstExonsOnlySizeSum),
        ("first_exons_only_unpredicted_size_sum", totalFirstExonsOnlySizeSum - successFirstExonsOnlySizeSum),
        ("last_exons_only_predicted_size_sum", successLastExonsOnlySizeSum),
        ("last_exons_only_unpredicted_size_sum", totalLastExonsOnlySizeSum - successLastExonsOnlySizeSum),
        ("single_exons_predicted_size_sum", successSingleExonsSizeSum),
        ("single_exons_unpredicted_size_sum", totalSingleExonsSizeSum - successSingleExonsSizeSum)
    ]
    strandedList = getStrandedVersionOfList(predictionDataGeneral,isForwardStrand)
    emptyList = getEmptyDf(not isForwardStrand,isForwardStrand,False,False)
    predictionData = predictionDataGeneral + strandedList + emptyList
    predictionDataDf = pd.DataFrame(predictionData, columns=["identifier", "value"])

    return predictionDataDf

def nucleotidesByGene(group):
    exons = group.loc[
        group["is_exon"], ["region_start", "region_end"]
    ].to_numpy(dtype=np.int64)

    if exons.size == 0:
        return ""

    starts = exons[:, 0]
    ends = exons[:, 1]

    lengths = ends - starts + 1
    totalLength = lengths.sum()
    nucleotides = np.empty(totalLength, dtype=np.int64)
    
    pos = 0
    for start, length in zip(starts, lengths):
        nucleotides[pos:pos+length] = np.arange(start, start+length)
        pos += length
    
    # Unique and sorted
    nucleotides = np.unique(nucleotides)
    
    return ";".join(map(str, nucleotides))

def getSingleColumnSet(df,mask,column):
    return set(df.loc[mask,column])

def getDoubleColumnSet(df,mask):
    setOfPairs = set(zip(
        df.loc[mask,'region_start'],
        df.loc[mask,'region_end']
    ))

    return setOfPairs

def singleColumnPredicted(df,mask,newColumn,comparisonColumn,comparisonList):
    df.loc[mask,newColumn] = df.loc[mask,comparisonColumn].isin(comparisonList)

def pairColumnPredicted(df,mask,listOfPairs):
    df.loc[mask, "region_predicted"] = (
        pd.MultiIndex.from_arrays([df.loc[mask, "region_start"], df.loc[mask, "region_end"]]).isin(listOfPairs)
    )

def setNumericalStatisticsSingleFile(baselineDf, candidateDf, chromosomePath, reader: Reader, groupType, isForwardStrand):
    mainDf, comparisonDf = [baselineDf, candidateDf] if groupType == "baseline" else [candidateDf, baselineDf]
    statisticChromosomePath = reader.getDefinedChromosomeStatistics(chromosomePath, groupType)
    statisticPath = reader.getFinalStatistics(groupType)
    nucleotidesPath = reader.getDefinedChromosomeNucleotides(chromosomePath, groupType)

    if mainDf is None:
        saveDf = getEmptyDf(True,True,True,True)
        nucleotideDf = pd.DataFrame(columns=["gene_id", "nucleotides"])
    else:
        df = mainDf.copy()
        df["donor_predicted"] = False
        df["acceptor_predicted"] = False
        df["region_predicted"] = False
        df["gene_predicted"] = False

        if comparisonDf is not None:

            intronMainPos, exonMainPos, startCodonMainPos, stopCodonMainPos = getMasks(df,['is_intron', 'is_exon', 'is_start_codon', 'is_stop_codon'])
            intronComparisionPos, exonComparisionPos, startCodonComparisionPos, stopCodonComparisionPos = getMasks(comparisonDf,['is_intron', 'is_exon', 'is_start_codon', 'is_stop_codon'])            

            comparisonIntronStart = getSingleColumnSet(comparisonDf,intronComparisionPos,'region_start')
            comparisonIntronEnd = getSingleColumnSet(comparisonDf,intronComparisionPos,'region_end')

            comparisonIntron = getDoubleColumnSet(comparisonDf,intronComparisionPos)
            comparisonExon = getDoubleColumnSet(comparisonDf,exonComparisionPos)
            comparisonStartCodon = getDoubleColumnSet(comparisonDf,startCodonComparisionPos)
            comparisonStopCodon = getDoubleColumnSet(comparisonDf,stopCodonComparisionPos)

            singleColumnPredicted(df,intronMainPos,'donor_predicted','region_start',comparisonIntronStart)
            singleColumnPredicted(df,intronMainPos,'acceptor_predicted','region_end',comparisonIntronEnd)

            pairColumnPredicted(df,intronMainPos,comparisonIntron)
            pairColumnPredicted(df,exonMainPos,comparisonExon)
            pairColumnPredicted(df,startCodonMainPos,comparisonStartCodon)
            pairColumnPredicted(df,stopCodonMainPos,comparisonStopCodon)

        nucleotideDf = (
            df
            .groupby("gene_id")
            .apply(nucleotidesByGene)
            .reset_index(name="nucleotides")
        )
        nucleotideDf["is_forward_strand"] = isForwardStrand
        saveDf = regionPredictionData(df,isForwardStrand)

    saveDf.to_csv(statisticChromosomePath, encoding='utf-8', index=False)
    saveDf.to_csv(statisticPath, encoding='utf-8', mode="a", index=False, header=False)    
    nucleotideDf.to_csv(nucleotidesPath, index=False)

def computeOverlapStats(df, refNucleotides):
    df["overlap"] = df["nucleotide_set"].apply(lambda s: bool(s & refNucleotides))
    total = len(df)
    not_ignored = df["overlap"].sum()
    ignored = total - not_ignored
    return total, not_ignored, ignored

def computeNucleotideSet(df):
    df["nucleotide_set"] = df["nucleotides"].apply(lambda x: set(map(int, x.split(";"))) if isinstance(x, str) else set())
    nucleotidesSet = set().union(*df["nucleotide_set"].tolist())

    return df, nucleotidesSet

def writeNucleotides(reader: Reader,chromosomePath,isForwardStrand):
    recallChromosomeCsv, precisionChromosomeCsv = reader.getDefinedChromosomeStatistics(chromosomePath, "baseline"), reader.getDefinedChromosomeStatistics(chromosomePath, "candidate")
    baselineNucleotidesCsv, candidateNucleotidesCsv = reader.getDefinedChromosomeNucleotides(chromosomePath, "baseline"), reader.getDefinedChromosomeNucleotides(chromosomePath, "candidate")
    recallCsv, precisionCsv = reader.getFinalStatistics("baseline"), reader.getFinalStatistics("candidate")
    
    baselineNucleotideDf = pd.read_csv(baselineNucleotidesCsv)
    candidateNucleotideDf = pd.read_csv(candidateNucleotidesCsv)


    baselineNucleotideDf, baselineNucleotides = computeNucleotideSet(baselineNucleotideDf)
    candidateNucleotideDf, candidateNucleotides = computeNucleotideSet(candidateNucleotideDf)

    recallGenesTotal, recallGenesNotIgnored, recallGenesIgnored = computeOverlapStats(baselineNucleotideDf,candidateNucleotides)
    precisionGenesTotal, precisionGenesNotIgnored, precisionGenesIgnored = computeOverlapStats(candidateNucleotideDf,baselineNucleotides)

    baselineNucleotidesSize = len(baselineNucleotides)
    candidateNucleotidesSize = len(candidateNucleotides)

    intersectionNucleotides = baselineNucleotides & candidateNucleotides
    intersectionNucleotidesSize = len(intersectionNucleotides)


    recallDataBase = [
        ("nucleotides_total", baselineNucleotidesSize),
        ("nucleotides_predicted", intersectionNucleotidesSize),
        ("nucleotides_unpredicted", (baselineNucleotidesSize - intersectionNucleotidesSize)),
        ("gene_total_size", recallGenesTotal),
        ("gene_ignored", recallGenesIgnored),
        ("gene_not_ignored", recallGenesNotIgnored)
    ]
    recallDataStranded = getStrandedVersionOfList(recallDataBase,isForwardStrand) 
    recallData = recallDataBase + recallDataStranded

    precisionDataBase = [
        ("nucleotides_total", candidateNucleotidesSize),
        ("nucleotides_predicted", intersectionNucleotidesSize),
        ("nucleotides_unpredicted", (candidateNucleotidesSize - intersectionNucleotidesSize)),
        ("gene_total_size", precisionGenesTotal),
        ("gene_ignored", precisionGenesIgnored),
        ("gene_not_ignored", precisionGenesNotIgnored)
    ]
    precisionDataStranded = getStrandedVersionOfList(precisionDataBase,isForwardStrand)
    precisionData = precisionDataBase + precisionDataStranded

    recallDataDf = pd.DataFrame(recallData, columns=reader.getStatisticsDf())
    precisionDataDf = pd.DataFrame(precisionData, columns=reader.getStatisticsDf())


    recallDataDf.to_csv(recallChromosomeCsv, encoding='utf-8', mode="a", index=False, header=False)
    recallDataDf.to_csv(recallCsv, encoding='utf-8', mode="a", index=False, header=False)

    precisionDataDf.to_csv(precisionChromosomeCsv, encoding='utf-8', mode="a", index=False, header=False)
    precisionDataDf.to_csv(precisionCsv, encoding='utf-8', mode="a", index=False, header=False)

def moveGeneTranscript(reader: Reader, chromosomePath):
    with open(reader.getDefinedChromosomeSingleGeneStringPath(chromosomePath),"r") as f_in, open(reader.getSingleGeneStringPath(),"a") as f_out:
        content = f_in.read()
        noHeaderContent = content.split("\n",1)[1]
        f_out.write(noHeaderContent)

def generateStatisticsPerFolder(reader: Reader, chromosomePath, isForwardStrand):
    baselineDf, candidateDf = reader.getProcessedData(chromosomePath, "baseline"), reader.getProcessedData(chromosomePath, "candidate")

    setNumericalStatisticsSingleFile(baselineDf, candidateDf, chromosomePath, reader, "baseline", isForwardStrand)
    setNumericalStatisticsSingleFile(baselineDf, candidateDf, chromosomePath, reader, "candidate", isForwardStrand)
    
    writeNucleotides(reader, chromosomePath, isForwardStrand)
    moveGeneTranscript(reader, chromosomePath)


def writeHeaders(reader: Reader):
    paths = [
        reader.getFinalStatistics("baseline"),
        reader.getFinalStatistics("candidate"),
        reader.getSingleGeneStringPath()
    ]
    headers = [
        reader.getStatisticsDf(),
        reader.getStatisticsDf(),
        reader.getGeneStringDfCols()
    ]

    for path, header in zip(paths, headers):
        with open(path, "w") as f:
            f.write(f"{','.join(header)}\n")

def getIntDivision(v1,v2,multiplier=1):
    if v2 == 0:
        return None
    
    return round(multiplier*(v1/v2),6)

def getGroupedDataDf(df,maxPos):
    identifiers = df["identifier"].unique()
    fullGrid = pd.MultiIndex.from_product([identifiers, range(maxPos + 1)], names=["identifier", "pos"]).to_frame(index=False)

    filledDf = fullGrid.merge(df, on=["identifier", "pos"], how="left")

    filledDf["value"] = filledDf["value"].fillna(0)

    return filledDf

def getTotals(df,cols):
    return [df.loc[df["identifier"] == col,"value"].iloc[0] for col in cols]

def calcAvg(dividend, divisor, multiplier=1):
    dividendSum = sum(dividend) if isinstance(dividend, list) else dividend
    divisorSum = sum(divisor) if isinstance(divisor, list) else divisor
    return round((multiplier*dividendSum) / divisorSum, 6) if divisorSum != 0 else None

def getDivisions(series,dictVar,listOfTriads,multiplier=1):
    for key, dividend, divisor in listOfTriads:
        par1Fwd = series[f"{dividend}_forward"]
        par2Fwd = series[f"{divisor}_forward"]
        fwdValue = calcAvg(par1Fwd,par2Fwd,multiplier)

        par1Rev = series[f"{dividend}_reverse"]
        par2Rev = series[f"{divisor}_reverse"]
        revValue = calcAvg(par1Rev,par2Rev,multiplier)

        par1Tot = series[f"{dividend}"]
        par2Tot = series[f"{divisor}"]
        totValue = calcAvg(par1Tot,par2Tot,multiplier)

        dictVar[f"{key}_forward"] = fwdValue
        dictVar[f"{key}_reverse"] = revValue
        dictVar[f"{key}"] = totValue

def getIntDivisions(dictVar,listOfTriads,multiplier=1):
    for key, dividend, divisor in listOfTriads:
        dictVar[key] = getIntDivision(dividend,divisor,multiplier)

def getSeries(reader: Reader, mainFilePath, isBaseline):
    #Genes Dfs
    genesDf = pd.read_csv(reader.getSingleGeneStringPath())
    genesDf = genesDf.loc[genesDf["is_baseline"] == isBaseline]

    sameStrandDf = genesDf.loc[~genesDf["gene_predicted"], ["chromosome_identifier", "gene_id", "is_forward_strand", "same_strand"]]
    sameStrandDf = sameStrandDf.groupby(["chromosome_identifier", "gene_id", "is_forward_strand"])["same_strand"].any().reset_index()
    sameStrandSeries = sameStrandDf.groupby(["is_forward_strand"]
    ).agg(
        gene_qtty=('same_strand', 'size'),
        correct_strand_sum=('same_strand', 'sum')
    )

    genesSizeSeries = genesDf.groupby(["gene_predicted","predicted","is_forward_strand"]
    ).agg(
        exon_qtty=('exon_qtty', 'size'),
        exon_sum=('exon_qtty', 'sum'),
        intron_retention_sum=('intron_retention_qtty', 'sum'),
        correct_frame_sum=('same_strand', 'sum')
    )

    genesGroupedGeneral = (
        genesDf.groupby(["chromosome_identifier","gene_id","is_forward_strand"])["gene_predicted"]
        .any()
    )

    genesPercentageSeries = genesGroupedGeneral.groupby("is_forward_strand").agg(
        predicted_sum='sum',
        gene_qtty='size'
    )

    genesPercentageSeries["unpredicted_sum"] = (genesPercentageSeries["gene_qtty"] - genesPercentageSeries["predicted_sum"])

    # Main series
    baseDf = pd.read_csv(mainFilePath)
    mainSeries = baseDf.groupby(["identifier"])["value"].sum()


    return sameStrandSeries, genesSizeSeries, genesPercentageSeries, mainSeries

def getSingleExonSumAndQtty(df, genePredicted, transcriptPredicted, isForwardStrand):
    key = (genePredicted, transcriptPredicted, isForwardStrand)

    if key not in df.index:
        return 0, 0, 0
    
    row = df.loc[key]
    return row["exon_qtty"], row["exon_sum"], row["intron_retention_sum"]

def getSameStrandOccurance(df, isForwardStrand):
    key = isForwardStrand

    if key not in df.index:
        return 0, 0

    row = df.loc[key]

    return row["correct_strand_sum"], row["gene_qtty"]

def getSingleGenePred(df, isForwardStrand):
    key = isForwardStrand

    if key not in df.index:
        return 0, 0, 0
    
    row = df.loc[key]

    return row["predicted_sum"], row["unpredicted_sum"], row["gene_qtty"]

def getExonSizeStatistics(df, generalStatistic):
    # Total number & quantity of exons in genes that were successfully predicted and transcripts that were successfully predicted
    fwdGenePredTranscPredExonQtty, fwdGenePredTranscPredExonSum, fwdGenePredTranscPredIntronRetentionSum  = getSingleExonSumAndQtty(df,True,True,True)
    revGenePredTranscPredExonQtty, revGenePredTranscPredExonSum, revGenePredTranscPredIntronRetentionSum  = getSingleExonSumAndQtty(df,True,True,False)
    # Exon Quantity
    fwdGenePredTranscPredExonAvg = calcAvg(fwdGenePredTranscPredExonSum,fwdGenePredTranscPredExonQtty)
    revGenePredTranscPredExonAvg = calcAvg(revGenePredTranscPredExonSum,revGenePredTranscPredExonQtty)
    totalGenePredTranscPredExonAvg = calcAvg([fwdGenePredTranscPredExonSum,revGenePredTranscPredExonSum],[fwdGenePredTranscPredExonQtty,revGenePredTranscPredExonQtty])
    # Exon of Intron Retention Quantity
    fwdGenePredTranscPredIntronRetentionAvg = calcAvg(fwdGenePredTranscPredIntronRetentionSum,fwdGenePredTranscPredExonQtty,100)
    revGenePredTranscPredIntronRetentionAvg = calcAvg(revGenePredTranscPredIntronRetentionSum,revGenePredTranscPredExonQtty,100)
    totalGenePredTranscPredIntronRetentionAvg = calcAvg([fwdGenePredTranscPredIntronRetentionSum,revGenePredTranscPredIntronRetentionSum],[fwdGenePredTranscPredExonQtty,revGenePredTranscPredExonQtty],100)

    # Total number & quantity of exons in genes that were successfully predicted and transcripts that were successfully predicted
    fwdGenePredTranscUnpredExonQtty, fwdGenePredTranscUnpredExonSum, fwdGenePredTranscUnpredIntronRetentionSum = getSingleExonSumAndQtty(df,True,False,True)
    revGenePredTranscUnpredExonQtty, revGenePredTranscUnpredExonSum, revGenePredTranscUnpredIntronRetentionSum = getSingleExonSumAndQtty(df,True,False,False)
    # Exon Quantity
    fwdGenePredTranscUnpredExonAvg = calcAvg(fwdGenePredTranscUnpredExonSum,fwdGenePredTranscUnpredExonQtty)
    revGenePredTranscUnpredExonAvg = calcAvg(revGenePredTranscUnpredExonSum,revGenePredTranscUnpredExonQtty)
    totalGenePredTranscUnpredExonAvg = calcAvg([fwdGenePredTranscUnpredExonSum,revGenePredTranscUnpredExonSum],[fwdGenePredTranscUnpredExonQtty,revGenePredTranscUnpredExonQtty])
    # Exon of Intron Retention Quantity
    fwdGenePredTranscUnpredIntronRetentionAvg = calcAvg(fwdGenePredTranscUnpredIntronRetentionSum,fwdGenePredTranscUnpredExonQtty,100)
    revGenePredTranscUnpredIntronRetentionAvg = calcAvg(revGenePredTranscUnpredIntronRetentionSum,revGenePredTranscUnpredExonQtty,100)
    totalGenePredTranscUnpredIntronRetentionAvg = calcAvg([fwdGenePredTranscUnpredIntronRetentionSum,revGenePredTranscUnpredIntronRetentionSum],[fwdGenePredTranscUnpredExonQtty,revGenePredTranscUnpredExonQtty],100)

    # Total number & quantity of exons in genes that were not predicted by any transcript
    fwdGeneUnpredTranscUnpredExonQtty, fwdGeneUnpredTranscUnpredExonSum, fwdGeneUnpredTranscUnpredIntronRetentionSum = getSingleExonSumAndQtty(df,False,False,True)
    revGeneUnpredTranscUnpredExonQtty, revGeneUnpredTranscUnpredExonSum, revGeneUnpredTranscUnpredIntronRetentionSum = getSingleExonSumAndQtty(df,False,False,False)
    # Exon Quantity
    fwdGeneUnpredTranscUnpredExonAvg = calcAvg(fwdGeneUnpredTranscUnpredExonSum,fwdGeneUnpredTranscUnpredExonQtty)
    revGeneUnpredTranscUnpredExonAvg = calcAvg(revGeneUnpredTranscUnpredExonSum,revGeneUnpredTranscUnpredExonQtty)
    totalGeneUnpredTranscUnpredExonAvg = calcAvg([fwdGeneUnpredTranscUnpredExonSum,revGeneUnpredTranscUnpredExonSum],[fwdGeneUnpredTranscUnpredExonQtty,revGeneUnpredTranscUnpredExonQtty])
    # Exon of Intron Retention Quantity
    fwdGeneUnpredTranscUnpredIntronRetentionAvg = calcAvg(fwdGeneUnpredTranscUnpredIntronRetentionSum,fwdGeneUnpredTranscUnpredExonQtty,100)
    revGeneUnpredTranscUnpredIntronRetentionAvg = calcAvg(revGeneUnpredTranscUnpredIntronRetentionSum,revGeneUnpredTranscUnpredExonQtty,100)
    totalGeneUnpredTranscUnpredIntronRetentionAvg = calcAvg([fwdGeneUnpredTranscUnpredIntronRetentionSum,revGeneUnpredTranscUnpredIntronRetentionSum],[fwdGeneUnpredTranscUnpredExonQtty,revGeneUnpredTranscUnpredExonQtty],100)


    #Updating generalStatistics
    # Exons: Predicted Genes and Predicted Transcripts
    generalStatistic["exons_in_predicted_genes_predicted_transcripts_forward"] = fwdGenePredTranscPredExonAvg
    generalStatistic["exons_in_predicted_genes_predicted_transcripts_reverse"] = revGenePredTranscPredExonAvg
    generalStatistic["exons_in_predicted_genes_predicted_transcripts"] = totalGenePredTranscPredExonAvg
    # Exons: Predicted Genes and Unpredicted Transcripts
    generalStatistic["exons_in_predicted_genes_unpredicted_transcripts_forward"] = fwdGenePredTranscUnpredExonAvg
    generalStatistic["exons_in_predicted_genes_unpredicted_transcripts_reverse"] = revGenePredTranscUnpredExonAvg
    generalStatistic["exons_in_predicted_genes_unpredicted_transcripts"] = totalGenePredTranscUnpredExonAvg
    # Exons: Unpredicted Genes and Unpredicted Transcripts
    generalStatistic["exons_in_unpredicted_genes_unpredicted_transcripts_forward"] = fwdGeneUnpredTranscUnpredExonAvg
    generalStatistic["exons_in_unpredicted_genes_unpredicted_transcripts_reverse"] = revGeneUnpredTranscUnpredExonAvg
    generalStatistic["exons_in_unpredicted_genes_unpredicted_transcripts"] = totalGeneUnpredTranscUnpredExonAvg
    # Intron Retention Exons: Predicted Genes and Predicted Transcripts
    generalStatistic["intron_retention_exons_in_predicted_genes_predicted_transcripts_forward_percentage"] = fwdGenePredTranscPredIntronRetentionAvg
    generalStatistic["intron_retention_exons_in_predicted_genes_predicted_transcripts_reverse_percentage"] = revGenePredTranscPredIntronRetentionAvg
    generalStatistic["intron_retention_exons_in_predicted_genes_predicted_transcripts_percentage"] = totalGenePredTranscPredIntronRetentionAvg
    # Intron Retention Exons: Predicted Genes and Unpredicted Transcripts
    generalStatistic["intron_retention_exons_in_predicted_genes_unpredicted_transcripts_forward_percentage"] = fwdGenePredTranscUnpredIntronRetentionAvg
    generalStatistic["intron_retention_exons_in_predicted_genes_unpredicted_transcripts_reverse_percentage"] = revGenePredTranscUnpredIntronRetentionAvg
    generalStatistic["intron_retention_exons_in_predicted_genes_unpredicted_transcripts_percentage"] = totalGenePredTranscUnpredIntronRetentionAvg
    # Intron Retention Exons: Unpredicted Genes and Unpredicted Transcripts
    generalStatistic["intron_retention_exons_in_unpredicted_genes_unpredicted_transcripts_forward_percentage"] = fwdGeneUnpredTranscUnpredIntronRetentionAvg
    generalStatistic["intron_retention_exons_in_unpredicted_genes_unpredicted_transcripts_reverse_percentage"] = revGeneUnpredTranscUnpredIntronRetentionAvg
    generalStatistic["intron_retention_exons_in_unpredicted_genes_unpredicted_transcripts_percentage"] = totalGeneUnpredTranscUnpredIntronRetentionAvg

    return

def getPredictionPercentageStatistics(df,generalStatistics):
    # Genes total size, total predicted, unpredicted and percentage

    # in the Forward Strand
    fwdGenePredSum, fwdGeneUnpredSum, fwdGeneQtty = getSingleGenePred(df,True) 
    fwdGenePredPerc = calcAvg(fwdGenePredSum,fwdGeneQtty,100)
    fwdGeneUnpredPerc = calcAvg(fwdGeneUnpredSum,fwdGeneQtty,100)

    # in the Reverse Strand
    revGenePredSum, revGeneUnpredSum, revGeneQtty = getSingleGenePred(df,False)
    revGenePredPerc = calcAvg(revGenePredSum,revGeneQtty,100)
    revGeneUnpredPerc = calcAvg(revGeneUnpredSum,revGeneQtty,100)

    #in Total
    totalGenePredPerc = calcAvg([fwdGenePredSum,revGenePredSum],[fwdGeneQtty,revGeneQtty],100)
    totalGeneUnpredPerc = calcAvg([fwdGeneUnpredSum,revGeneUnpredSum],[fwdGeneQtty,revGeneQtty],100)

    #Updating generalStatistics
    # Genes Predicted
    generalStatistics["genes_predicted_percentage_forward"] = fwdGenePredPerc
    generalStatistics["genes_predicted_percentage_reverse"] = revGenePredPerc
    generalStatistics["genes_predicted_percentage"] = totalGenePredPerc
    # Genes Unpredicted
    generalStatistics["genes_unpredicted_percentage_forward"] = fwdGeneUnpredPerc
    generalStatistics["genes_unpredicted_percentage_reverse"] = revGeneUnpredPerc
    generalStatistics["genes_unpredicted_percentage"] = totalGeneUnpredPerc

    return

def getSameStrandStatistics(df, generalStatistics, isBaseline):
    fwdSameStrandCount, fwdGeneQtty = getSameStrandOccurance(df, True)
    revSameStrandCount, revGeneQtty = getSameStrandOccurance(df, False)

    fwdUnpredGeneSameStrandPerc = calcAvg(fwdSameStrandCount, fwdGeneQtty,100)
    revUnpredGeneSameStrandPerc = calcAvg(revSameStrandCount, revGeneQtty,100)
    totalUnpredGeneSameStrandPerc = calcAvg([fwdSameStrandCount, revSameStrandCount], [fwdGeneQtty, revGeneQtty],100)

    generalStatistics["unpredicted_genes_with_same_strand_percentage_forward"] = fwdUnpredGeneSameStrandPerc if isBaseline else 0.0
    generalStatistics["unpredicted_genes_with_same_strand_percentage_reverse"] = revUnpredGeneSameStrandPerc if isBaseline else 0.0
    generalStatistics["unpredicted_genes_with_same_strand_percentage"] = totalUnpredGeneSameStrandPerc if isBaseline else 0.0

    return

def generateGeneralStatistics(reader: Reader, filePath, statPath, isBaseline):
    generalStatistic = {}
    sameStrandSeries, genesSizeSeries, genesPercentageSeries, mainSeries = getSeries(reader, filePath, isBaseline)

    getSameStrandStatistics(sameStrandSeries, generalStatistic, isBaseline)
    getExonSizeStatistics(genesSizeSeries, generalStatistic)
    getPredictionPercentageStatistics(genesPercentageSeries, generalStatistic)

    getDivisions(
        mainSeries,
        generalStatistic,
        [
            ["genes_ignored_percentage","gene_ignored","gene_total_size"],
            ["genes_not_ignored_percentage","gene_not_ignored","gene_total_size"],
            ["nucleotides_predicted_percentage","nucleotides_predicted","nucleotides_total"],
            ["nucleotides_unpredicted_percentage","nucleotides_unpredicted","nucleotides_total"],
            ["exons_predicted_percentage","exons_predicted","exons_total"],
            ["exons_unpredicted_percentage","exons_unpredicted","exons_total"],
            ["introns_predicted_percentage","introns_predicted","introns_total"],
            ["introns_unpredicted_percentage","introns_unpredicted","introns_total"],
            ["start_codons_predicted_percentage","start_codons_predicted","start_codons_total"],
            ["start_codons_unpredicted_percentage","start_codons_unpredicted","start_codons_total"],
            ["stop_codons_predicted_percentage","stop_codons_predicted","stop_codons_total"],
            ["stop_codons_unpredicted_percentage","stop_codons_unpredicted","stop_codons_total"],
            ["first_exons_only_predicted_percentage","first_exons_only_predicted","first_exons_only_total"],
            ["first_exons_only_unpredicted_percentage","first_exons_only_unpredicted","first_exons_only_total"],
            ["last_exons_only_predicted_percentage","last_exons_only_predicted","last_exons_only_total"],
            ["last_exons_only_unpredicted_percentage","last_exons_only_unpredicted","last_exons_only_total"],
            ["single_exons_predicted_percentage","single_exons_predicted","single_exons_total"],
            ["single_exons_unpredicted_percentage","single_exons_unpredicted","single_exons_total"],
            ["donors_predicted_percentage","donors_predicted","donors_total"],
            ["donors_unpredicted_percentage","donors_unpredicted","donors_total"],
            ["acceptors_predicted_percentage","acceptors_predicted","acceptors_total"],
            ["acceptors_unpredicted_percentage","acceptors_unpredicted","acceptors_total"]
        ],
        100
    )

    getDivisions(
        mainSeries,
        generalStatistic,
        [
            ["exons_total_size_avg","exons_total_size_sum","exons_total"],
            ["exons_predicted_size_avg","exons_predicted_size_sum","exons_predicted"],
            ["exons_unpredicted_size_avg","exons_unpredicted_size_sum","exons_unpredicted"],
            ["introns_total_size_avg","introns_total_size_sum","introns_total"],
            ["introns_predicted_size_avg","introns_predicted_size_sum","introns_predicted"],
            ["introns_unpredicted_size_avg","introns_unpredicted_size_sum","introns_unpredicted"],
            ["first_exons_only_total_size_avg","first_exons_only_total_size_sum","first_exons_only_total"],
            ["first_exons_only_predicted_size_avg","first_exons_only_predicted_size_sum","first_exons_only_predicted"],
            ["first_exons_only_unpredicted_size_avg","first_exons_only_unpredicted_size_sum","first_exons_only_unpredicted"],
            ["last_exons_only_total_size_avg","last_exons_only_total_size_sum","last_exons_only_total"],
            ["last_exons_only_predicted_size_avg","last_exons_only_predicted_size_sum","last_exons_only_predicted"],
            ["last_exons_only_unpredicted_size_avg","last_exons_only_unpredicted_size_sum","last_exons_only_unpredicted"],
            ["single_exons_total_size_avg","single_exons_total_size_sum","single_exons_total"],
            ["single_exons_predicted_size_avg","single_exons_predicted_size_sum","single_exons_predicted"],
            ["single_exons_unpredicted_size_avg","single_exons_unpredicted_size_sum","single_exons_unpredicted"]
        ]
    )

    with open(statPath,"w") as outputFile:
        json.dump(generalStatistic, outputFile, indent=4)

def generateMultipleGeneralStatistics(reader: Reader):
    filePathRecall, filePathPrecision = reader.getFinalStatistics("baseline"), reader.getFinalStatistics("candidate")
    statPathRecall, statPathPrecision = reader.getFinalStatisticsJson("baseline"), reader.getFinalStatisticsJson("candidate")

    generateGeneralStatistics(reader, filePathRecall, statPathRecall, True)
    generateGeneralStatistics(reader, filePathPrecision, statPathPrecision, False)


def statisticalAnalysis(reader: Reader):
    writeHeaders(reader)

    chromosomeFolders = reader.getChromosomeFoldersList()
    for chromosomePath in chromosomeFolders:
        isForwardStrand = chromosomePath.split("__")[-1] == "forward_strand"
        generateStatisticsPerFolder(reader, chromosomePath, isForwardStrand)
    generateMultipleGeneralStatistics(reader)
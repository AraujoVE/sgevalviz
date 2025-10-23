import pandas as pd
import json
import os
from sgevalviz.utils import *
import itertools
import numpy as np

def getChromosomeFolders(saveFilesBasePath):
    folder = f"{saveFilesBasePath}chromosomeCSVs"
    subfolders = [f.path for f in os.scandir(folder) if f.is_dir()]
    
    return subfolders

def getMasks(df,maskStrList):
    return [df[col].to_numpy() for col in maskStrList]

def groupDf(df,groupbyList,aggTerm):
    return df.groupby(groupbyList)[aggTerm].any().reset_index()

def addSizeParams(df,regionSize):
    df["region_length"] = (df["region_end"] - df["region_start"]).abs() + 1
    df["size_group"] = df['region_length'] // regionSize
    df["is_long"] = df["region_length"] > 300

def intersectMasks(mainMask, secondaryMaskList):
    return [(mainMask & secondaryMask) for secondaryMask in secondaryMaskList]

def countTrueValues(masks):
    return [mask.sum() for mask in masks]

def countTotalSizeSum(regionLength, maskList):
    return [(regionLength[mask]).sum() for mask in maskList]

def regionPredictionData(mainDf,regionSize):
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
    addSizeParams(mainDfDedup,regionSize)
    addSizeParams(firstsDf,regionSize)

    # Region Length
    mainDfDedupRegionLength = mainDfDedup["region_length"].to_numpy()
    firstsDfRegionLength = firstsDf["region_length"].to_numpy()

    # Base Masks
    firstDfPredicted, isFirstExon, isLastExon = getMasks(firstsDf,["region_predicted","is_first_exon","is_last_exon"]) 
    mainRegionPredicted, isExon, isIntron, isLong, isStartCodon, isStopCodon = getMasks(mainDfDedup,["region_predicted","is_exon","is_intron","is_long","is_start_codon","is_stop_codon"])
    isDonorPredicted = getMasks(donorDf,["donor_predicted"])[0]
    isAcceptorPredicted = getMasks(acceptorDf,["acceptor_predicted"])[0]

    # Secondary Masks
    isShort = ~isLong
    isShortExon, isLongExon = intersectMasks(isExon, [isShort,isLong])
    isShortIntron, isLongIntron = intersectMasks(isIntron, [isShort,isLong])
    isExonPredicted, isIntronPredicted, isShortExonPredicted, isLongExonPredicted, isShortIntronPredicted, isLongIntronPredicted, isStartCodonPredicted, isStopCodonPredicted = intersectMasks(mainRegionPredicted,[isExon,isIntron,isShortExon,isLongExon,isShortIntron,isLongIntron,isStartCodon,isStopCodon])
    isFirstExonPredicted, isLastExonPredicted = intersectMasks(firstDfPredicted,[isFirstExon,isLastExon])

    # Count Quantities
    totalDonors, totalAcceptors = len(donorDf), len(acceptorDf)
    totalExons, totalShortExons, totalLongExons, totalIntrons, totalShortIntrons, totalLongIntrons, totalStartCodons, totalStopCodons, totalFirstExons, totalLastExons, successExons, successShortExons, successLongExons, successIntrons, successShortIntrons, successLongIntrons, successStartCodons, successStopCodons, successFirstExons, successLastExons, successDonors, successAcceptors = countTrueValues([isExon,isShortExon,isLongExon,isIntron,isShortIntron,isLongIntron,isStartCodon,isStopCodon,isFirstExon,isLastExon,isExonPredicted,isShortExonPredicted,isLongExonPredicted,isIntronPredicted,isShortIntronPredicted,isLongIntronPredicted,isStartCodonPredicted,isStopCodonPredicted,isFirstExonPredicted,isLastExonPredicted,isDonorPredicted,isAcceptorPredicted])

    # Count Sizes
    totalExonsSizeSum, totalShortExonsSizeSum, totalLongExonsSizeSum, totalIntronsSizeSum, totalShortIntronsSizeSum, totalLongIntronsSizeSum, successExonsSizeSum, successShortExonsSizeSum, successLongExonsSizeSum, successIntronsSizeSum, successShortIntronsSizeSum, successLongIntronsSizeSum = countTotalSizeSum(mainDfDedupRegionLength,[isExon,isShortExon,isLongExon,isIntron,isShortIntron,isLongIntron,isExonPredicted,isShortExonPredicted,isLongExonPredicted,isIntronPredicted,isShortIntronPredicted,isLongIntronPredicted])
    totalFirstExonsSizeSum, totalLastExonsSizeSum, successFirstExonsSizeSum, successLastExonsSizeSum = countTotalSizeSum(firstsDfRegionLength,[isFirstExon,isLastExon,isFirstExonPredicted,isLastExonPredicted])

    # Write Data
    predictionData = [
        ("exons_total", -1, totalExons),
        ("short_exons_total", -1, totalShortExons),
        ("long_exons_total", -1, totalLongExons),
        ("introns_total", -1, totalIntrons),
        ("short_introns_total", -1, totalShortIntrons),
        ("long_introns_total", -1, totalLongIntrons),
        ("start_codons_total", -1, totalStartCodons),
        ("stop_codons_total", -1, totalStopCodons),
        ("first_exons_total", -1, totalFirstExons),
        ("last_exons_total", -1, totalLastExons),
        ("donors_total", -1, totalDonors),
        ("acceptors_total", -1, totalAcceptors),
        ("exons_predicted", -1, successExons),
        ("exons_unpredicted", -1, totalExons - successExons),
        ("short_exons_predicted", -1, successShortExons),
        ("short_exons_unpredicted", -1, totalShortExons - successShortExons),
        ("long_exons_predicted", -1, successLongExons),
        ("long_exons_unpredicted", -1, totalLongExons - successLongExons),
        ("introns_predicted", -1, successIntrons),
        ("introns_unpredicted", -1, totalIntrons - successIntrons),
        ("short_introns_predicted", -1, successShortIntrons),
        ("short_introns_unpredicted", -1, totalShortIntrons - successShortIntrons),
        ("long_introns_predicted", -1, successLongIntrons),
        ("long_introns_unpredicted", -1, totalLongIntrons - successLongIntrons),
        ("start_codons_predicted", -1, successStartCodons),
        ("start_codons_unpredicted", -1, totalStartCodons - successStartCodons),
        ("stop_codons_predicted", -1, successStopCodons),
        ("stop_codons_unpredicted", -1, totalStopCodons - successStopCodons),
        ("first_exons_predicted", -1, successFirstExons),
        ("first_exons_unpredicted", -1, totalFirstExons - successFirstExons),
        ("last_exons_predicted", -1, successLastExons),
        ("last_exons_unpredicted", -1, totalLastExons - successLastExons),
        ("donors_predicted", -1, successDonors),
        ("donors_unpredicted", -1, totalDonors - successDonors),
        ("acceptors_predicted", -1, successAcceptors),
        ("acceptors_unpredicted", -1, totalAcceptors - successAcceptors),
        ("exons_total_size_sum", -1, totalExonsSizeSum),
        ("short_exons_total_size_sum", -1, totalShortExonsSizeSum),
        ("long_exons_total_size_sum", -1, totalLongExonsSizeSum),
        ("introns_total_size_sum", -1, totalIntronsSizeSum),
        ("short_introns_total_size_sum", -1, totalShortIntronsSizeSum),
        ("long_introns_total_size_sum", -1, totalLongIntronsSizeSum),
        ("first_exons_total_size_sum", -1, totalFirstExonsSizeSum),
        ("last_exons_total_size_sum", -1, totalLastExonsSizeSum),
        ("exons_predicted_size_sum", -1, successExonsSizeSum),
        ("exons_unpredicted_size_sum", -1, totalExonsSizeSum - successExonsSizeSum),
        ("short_exons_predicted_size_sum", -1, successShortExonsSizeSum),
        ("short_exons_unpredicted_size_sum", -1, totalShortExonsSizeSum - successShortExonsSizeSum),
        ("long_exons_predicted_size_sum", -1, successLongExonsSizeSum),
        ("long_exons_unpredicted_size_sum", -1, totalLongExonsSizeSum - successLongExonsSizeSum),
        ("introns_predicted_size_sum", -1, successIntronsSizeSum),
        ("introns_unpredicted_size_sum", -1, totalIntronsSizeSum - successIntronsSizeSum),
        ("short_introns_predicted_size_sum", -1, successShortIntronsSizeSum),
        ("short_introns_unpredicted_size_sum", -1, totalShortIntronsSizeSum - successShortIntronsSizeSum),
        ("long_introns_predicted_size_sum", -1, successLongIntronsSizeSum),
        ("long_introns_unpredicted_size_sum", -1, totalLongIntronsSizeSum - successLongIntronsSizeSum),
        ("first_exons_predicted_size_sum", -1, successFirstExonsSizeSum),
        ("first_exons_unpredicted_size_sum", -1, totalFirstExonsSizeSum - successFirstExonsSizeSum),
        ("last_exons_predicted_size_sum", -1, successLastExonsSizeSum),
        ("last_exons_unpredicted_size_sum", -1, totalLastExonsSizeSum - successLastExonsSizeSum)
    ]

    maxGroupSize = mainDfDedup["size_group"].max()

    for i in range(maxGroupSize+1):
        isFromGroupMain = (mainDfDedup["size_group"] == i).to_numpy(dtype=bool)
        isFromGroupFirstLast = (firstsDf["size_group"] == i).to_numpy(dtype=bool)

        # Masks
        isExonGroup, isIntronGroup, isExonPredictedGroup, isIntronPredictedGroup = intersectMasks(isFromGroupMain,[isExon,isIntron,isExonPredicted,isIntronPredicted])
        isFirstExonGroup, isLastExonGroup, isFirstExonPredictedGroup, isLastExonPredictedGroup = intersectMasks(isFromGroupFirstLast,[isFirstExon,isLastExon,isFirstExonPredicted,isLastExonPredicted])

        # Count Quantities
        totalExons, totalIntrons, totalExonsPredicted, totalIntronsPredicted, totalFirstExons, totalLastExons, totalFirstExonsPredicted, totalLastExonsPredicted = countTrueValues([isExonGroup,isIntronGroup,isExonPredictedGroup,isIntronPredictedGroup,isFirstExonGroup,isLastExonGroup,isFirstExonPredictedGroup,isLastExonPredictedGroup])


        # Count Sizes
        totalExonsSizeSum, totalIntronsSizeSum, totalExonsPredictedSizeSum, totalIntronsPredictedSizeSum = countTotalSizeSum(mainDfDedupRegionLength,[isExonGroup,isIntronGroup,isExonPredictedGroup,isIntronPredictedGroup])
        totalFirstExonsSizeSum, totalLastExonsSizeSum, totalFirstExonsPredictedSizeSum, totalLastExonsPredictedSizeSum = countTotalSizeSum(firstsDfRegionLength,[isFirstExonGroup,isLastExonGroup,isFirstExonPredictedGroup,isLastExonPredictedGroup])


        # Write Data
        partialPredictionData = [
            ("exons_total", i, totalExons),
            ("introns_total", i, totalIntrons),
            ("first_exons_total", i, totalFirstExons),
            ("last_exons_total", i, totalLastExons),
            ("exons_predicted", i, totalExonsPredicted),
            ("exons_unpredicted", i, totalExons - totalExonsPredicted),
            ("introns_predicted", i, totalIntronsPredicted),
            ("introns_unpredicted", i, totalIntrons - totalIntronsPredicted),
            ("first_exons_predicted", i, totalFirstExonsPredicted),
            ("first_exons_unpredicted", i, totalFirstExons - totalFirstExonsPredicted),
            ("last_exons_predicted", i, totalLastExonsPredicted),
            ("last_exons_unpredicted", i, totalLastExons - totalLastExonsPredicted),
            ("exons_total_size_sum", i, totalExonsSizeSum),
            ("introns_total_size_sum", i, totalIntronsSizeSum),
            ("first_exons_total_size_sum", i, totalFirstExonsSizeSum),
            ("last_exons_total_size_sum", i, totalLastExonsSizeSum),
            ("exons_predicted_size_sum", i, totalExonsPredictedSizeSum),
            ("exons_unpredicted_size_sum", i, totalExonsSizeSum - totalExonsPredictedSizeSum),
            ("introns_predicted_size_sum", i, totalIntronsPredictedSizeSum),
            ("introns_unpredicted_size_sum", i, totalIntronsSizeSum - totalIntronsPredictedSizeSum),
            ("first_exons_predicted_size_sum", i, totalFirstExonsPredictedSizeSum),
            ("first_exons_unpredicted_size_sum", i, totalFirstExonsSizeSum - totalFirstExonsPredictedSizeSum),
            ("last_exons_predicted_size_sum", i, totalLastExonsPredictedSizeSum),
            ("last_exons_unpredicted_size_sum", i, totalLastExonsSizeSum - totalLastExonsPredictedSizeSum)
        ]

        predictionData += partialPredictionData
    
    predictionDataDf = pd.DataFrame(predictionData, columns=["identifier", "pos", "value"])

    return predictionDataDf

def emptyDf():
    predictionData = [(par,-1,0) for par in [
        "exons_total",
        "short_exons_total",
        "long_exons_total",
        "introns_total",
        "short_introns_total",
        "long_introns_total",
        "start_codons_total",
        "stop_codons_total",
        "first_exons_total",
        "last_exons_total",
        "donors",
        "acceptors",
        "exons_predicted",
        "exons_unpredicted",
        "short_exons_predicted",
        "short_exons_unpredicted",
        "long_exons_predicted",
        "long_exons_unpredicted",
        "introns_predicted",
        "introns_unpredicted",
        "short_introns_predicted",
        "short_introns_unpredicted",
        "long_introns_predicted",
        "long_introns_unpredicted",
        "start_codons_predicted",
        "start_codons_unpredicted",
        "stop_codons_predicted",
        "stop_codons_unpredicted",
        "first_exons_predicted",
        "first_exons_unpredicted",
        "last_exons_predicted",
        "last_exons_unpredicted",
        "donors_predicted",
        "donors_unpredicted",
        "acceptors_predicted",
        "acceptors_unpredicted",
        "exons_total_size_sum",
        "short_exons_total_size_sum",
        "long_exons_total_size_sum",
        "introns_total_size_sum",
        "short_introns_total_size_sum",
        "long_introns_total_size_sum",
        "first_exons_total_size_sum",
        "last_exons_total_size_sum",
        "exons_predicted_size_sum",
        "exons_unpredicted_size_sum",
        "short_exons_predicted_size_sum",
        "short_exons_unpredicted_size_sum",
        "long_exons_predicted_size_sum",
        "long_exons_unpredicted_size_sum",
        "introns_predicted_size_sum",
        "introns_unpredicted_size_sum",
        "short_introns_predicted_size_sum",
        "short_introns_unpredicted_size_sum",
        "long_introns_predicted_size_sum",
        "long_introns_unpredicted_size_sum",
        "first_exons_predicted_size_sum",
        "first_exons_unpredicted_size_sum",
        "last_exons_predicted_size_sum",
        "last_exons_unpredicted_size_sum"
    ]]

    predictionDataDf = pd.DataFrame(predictionData, columns=["identifier", "pos", "value"])

    return predictionDataDf

def nucleotidesByGene(group):
    exons = group.loc[group["is_exon"], ["region_start", "region_end"]].to_numpy()
    
    if exons.size == 0:
        return ""
    
    starts = exons[:, 0]
    ends = exons[:, 1]
    
    # Compute lengths of each exon
    lengths = ends - starts + 1
    totalLength = lengths.sum()
    
    # Preallocate array for all nucleotide positions
    nucleotides = np.empty(totalLength, dtype=int)
    
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

def setNumericalStatisticsSingleFile(mainDf,comparisonDf,saveFilePath,saveFilePathFinal,regionSize,nucleotidesFilePath):

    if mainDf is None:
        saveDf = emptyDf()
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
        saveDf = regionPredictionData(df,regionSize)

    saveDf.to_csv(saveFilePath, encoding='utf-8', index=False)
    saveDf.to_csv(saveFilePathFinal, encoding='utf-8', mode="a", index=False, header=False)    
    nucleotideDf.to_csv(nucleotidesFilePath, index=False)

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

def writeNucleotides(saveFilesBasePath,chromosomeFolder):
    recallChromosomeCsv = f"{chromosomeFolder}/recallStatistics.csv"
    recallNucleotidesCsv = f"{chromosomeFolder}/recallNucleotides.csv"
    recallCsv = f"{saveFilesBasePath}chromosomeCSVs/recallStatistics.csv"

    precisionChromosomeCsv = f"{chromosomeFolder}/precisionStatistics.csv"
    precisionNucleotidesCsv = f"{chromosomeFolder}/precisionNucleotides.csv"
    precisionCsv = f"{saveFilesBasePath}chromosomeCSVs/precisionStatistics.csv"

    
    recallNucleotideDf = pd.read_csv(recallNucleotidesCsv)
    precisionNucleotideDf = pd.read_csv(precisionNucleotidesCsv)


    recallNucleotideDf, recallNucleotides = computeNucleotideSet(recallNucleotideDf)
    precisionNucleotideDf, precisionNucleotides = computeNucleotideSet(precisionNucleotideDf)

    recallGenesTotal, recallGenesNotIgnored, recallGenesIgnored = computeOverlapStats(recallNucleotideDf,precisionNucleotides)
    precisionGenesTotal, precisionGenesNotIgnored, precisionGenesIgnored = computeOverlapStats(precisionNucleotideDf,recallNucleotides)

    recallNucleotidesSize = len(recallNucleotides)
    precisionNucleotidesSize = len(precisionNucleotides)

    intersectionNucleotides = recallNucleotides & precisionNucleotides
    intersectionNucleotidesSize = len(intersectionNucleotides)


    recallData = [
        ("nucleotides_total", -1, recallNucleotidesSize),
        ("nucleotides_predicted", -1, intersectionNucleotidesSize),
        ("nucleotides_unpredicted", -1, (recallNucleotidesSize - intersectionNucleotidesSize)),
        ("gene_total_size", -1, recallGenesTotal),
        ("gene_ignored", -1, recallGenesIgnored),
        ("gene_not_ignored", -1, recallGenesNotIgnored)
    ]

    precisionData = [
        ("nucleotides_total", -1, precisionNucleotidesSize),
        ("nucleotides_predicted", -1, intersectionNucleotidesSize),
        ("nucleotides_unpredicted", -1, (precisionNucleotidesSize - intersectionNucleotidesSize)),
        ("gene_total_size", -1, precisionGenesTotal),
        ("gene_ignored", -1, precisionGenesIgnored),
        ("gene_not_ignored", -1, precisionGenesNotIgnored)
    ]

    recallDataDf = pd.DataFrame(recallData, columns=["identifier", "pos", "value"])
    precisionDataDf = pd.DataFrame(precisionData, columns=["identifier", "pos", "value"])


    recallDataDf.to_csv(recallChromosomeCsv, encoding='utf-8', mode="a", index=False, header=False)
    recallDataDf.to_csv(recallCsv, encoding='utf-8', mode="a", index=False, header=False)

    precisionDataDf.to_csv(precisionChromosomeCsv, encoding='utf-8', mode="a", index=False, header=False)
    precisionDataDf.to_csv(precisionCsv, encoding='utf-8', mode="a", index=False, header=False)

def generateStatisticsPerFolder(saveFilesBasePath,sf,regionSize):
    hasBaseline = os.path.exists(f"{sf}/processedBaselineFile.csv")
    hasCandidate = os.path.exists(f"{sf}/processedCandidateFile.csv")

    baselineDf = pd.read_csv(f"{sf}/processedBaselineFile.csv") if hasBaseline else None
    candidateDf = pd.read_csv(f"{sf}/processedCandidateFile.csv") if hasCandidate else None

    setNumericalStatisticsSingleFile(baselineDf,candidateDf,f"{sf}/recallStatistics.csv",f"{saveFilesBasePath}chromosomeCSVs/recallStatistics.csv",regionSize,f"{sf}/recallNucleotides.csv")
    setNumericalStatisticsSingleFile(candidateDf,baselineDf,f"{sf}/precisionStatistics.csv",f"{saveFilesBasePath}chromosomeCSVs/precisionStatistics.csv",regionSize,f"{sf}/precisionNucleotides.csv")    
    
    writeNucleotides(saveFilesBasePath,sf)

def writeHeaders(saveFilesBasePath,paths,contents):
    for i in range(len(paths)):
        with open(f"{saveFilesBasePath}{paths[i]}",'w') as f:
            f.write(f"{contents[i]}\n")

def moveGeneTranscript(saveFilesBasePath,sf):
    with open(f"{sf}/gene_transcript_predicted.csv","r") as f_in, open(f"{saveFilesBasePath}chromosomeCSVs/gene_transcript_predicted.csv","a") as f_out:
        content = f_in.read()
        noHeaderContent = content.split("\n",1)[1]
        f_out.write(noHeaderContent)

def getIntDivision(v1,v2,multiplier=1):
    if v2 == 0:
        return None
    
    return round(multiplier*(v1/v2),6)

def getTotal(df,p):
    return df.loc[df["identifier"] == p,"value"].sum()

def getSizeStatistic(df, isGenePredicted, isTranscriptPredicted):
    matches = df.loc[
        (df["gene_predicted"] == isGenePredicted) &
        (df["predicted"] == isTranscriptPredicted),
        "exon_qtty"
    ]
    
    if matches.empty:
        return None  # or 0, or float("nan"), depending on what makes sense
    
    return matches.iat[0]

def getGroupedDataDf(df,maxPos):
    identifiers = df["identifier"].unique()
    fullGrid = pd.MultiIndex.from_product([identifiers, range(maxPos + 1)], names=["identifier", "pos"]).to_frame(index=False)

    filledDf = fullGrid.merge(df, on=["identifier", "pos"], how="left")

    filledDf["value"] = filledDf["value"].fillna(0)

    return filledDf

def breakDf(filePath):
    baseDf = pd.read_csv(filePath)
    baseDf = baseDf.groupby(["identifier","pos"], as_index=False)["value"].sum()
    baseGroupedDf = baseDf.loc[baseDf["pos"] != -1].copy()
    ungroupedDf = baseDf.loc[baseDf["pos"] == -1].copy()
    ungroupedDf.drop(columns=["pos"], inplace=True)
    maxPos = baseGroupedDf["pos"].max()
    groupedDf = getGroupedDataDf(baseGroupedDf,maxPos)

    return groupedDf, ungroupedDf, maxPos

def getTotals(df,cols):
    return [getTotal(df,col) for col in cols]

def getDivisions(df,dictVar,listOfTriads,multiplier=1):
    valueMap = df.set_index("identifier")["value"]
    for key, dividend, divisor in listOfTriads:
        par1 = valueMap[dividend]
        par2 = valueMap[divisor]
        dictVar[key] = None if par2 == 0 else round(multiplier * (par1 / par2), 6)

def getIntDivisions(dictVar,listOfTriads,multiplier=1):
    for key, dividend, divisor in listOfTriads:
        dictVar[key] = getIntDivision(dividend,divisor,multiplier)

def getSizeStatistics(df,dictVar,listOfTriads):
    for key, isGenePredicted, isTranscriptPredicted in listOfTriads:
        dictVar[key] = getSizeStatistic(df,isGenePredicted,isTranscriptPredicted)

def initializeStatisticGroup(dictVar,groupIter,regionSize):
    dictVar["groups"][str(groupIter)] = {"size_range": f"{groupIter*regionSize} - {groupIter*regionSize + regionSize - 1}", "data": {}}

def getGenesData(saveFilesBasePath,isBaseline):
    genesDf = pd.read_csv(f"{saveFilesBasePath}chromosomeCSVs/gene_transcript_predicted.csv")
    genesDf = genesDf.loc[genesDf["is_baseline"] == isBaseline]
    
    genesSizeDf = (genesDf.groupby(["gene_predicted","predicted"])["exon_qtty"].mean().reset_index())
    genesGrouped = (genesDf.groupby(["chromosome_identifier","gene_id"])["gene_predicted"].any().reset_index())

    totalGenes = len(genesGrouped)

    predictedPercentage, unpredictedPercentage = None, None

    if totalGenes != 0:
        predictedGenes = genesGrouped["gene_predicted"].sum()
        predictedPercentage = round(100*(predictedGenes / totalGenes),6)
        unpredictedPercentage = round(100 - predictedPercentage,6)

    return genesSizeDf, predictedPercentage, unpredictedPercentage    

def generateGeneralStatistics(saveFilesBasePath,filePath,statPath,regionSize,isBaseline):
    generalStatistic = {}

    genesSizeDf, predictedPercentage, unpredictedPercentage = getGenesData(saveFilesBasePath,isBaseline)

    groupedDf, ungroupedDf, maxPos = breakDf(filePath)

    exonsQtty, exonsPredictedQtty, exonsUnpredictedQtty, firstExonsQtty, firstExonsPredictedQtty, firstExonsUnpredictedQtty, lastExonsQtty, lastExonsPredictedQtty, lastExonsUnpredictedQtty, intronsQtty, intronsPredictedQtty, intronsUnpredictedQtty = getTotals(ungroupedDf,["exons_total","exons_predicted","exons_unpredicted","first_exons_total","first_exons_predicted","first_exons_unpredicted","last_exons_total","last_exons_predicted","last_exons_unpredicted","introns_total","introns_predicted","introns_unpredicted"])


    getSizeStatistics(
        genesSizeDf,
        generalStatistic,
        [
            ["exons_in_predicted_genes_predicted_transcripts",True,True],
            ["exons_in_predicted_genes_unpredicted_transcripts",True,False],
            ["exons_in_unpredicted_genes",False,False]
        ]
    )

    generalStatistic["genes_predicted_percentage"] = predictedPercentage
    generalStatistic["genes_unpredicted_percentage"] = unpredictedPercentage

    getDivisions(
        ungroupedDf,
        generalStatistic,
        [
            ["genes_ignored_percentage","gene_ignored","gene_total_size"],
            ["genes_not_ignored_percentage","gene_not_ignored","gene_total_size"],
            ["nucleotides_predicted_percentage","nucleotides_predicted","nucleotides_total"],
            ["nucleotides_unpredicted_percentage","nucleotides_unpredicted","nucleotides_total"],
            ["exons_predicted_percentage","exons_predicted","exons_total"],
            ["exons_unpredicted_percentage","exons_unpredicted","exons_total"],
            ["short_exons_predicted_percentage","short_exons_predicted","short_exons_total"],
            ["short_exons_unpredicted_percentage","short_exons_unpredicted","short_exons_total"],
            ["long_exons_predicted_percentage","long_exons_predicted","long_exons_total"],
            ["long_exons_unpredicted_percentage","long_exons_unpredicted","long_exons_total"],
            ["introns_predicted_percentage","introns_predicted","introns_total"],
            ["introns_unpredicted_percentage","introns_unpredicted","introns_total"],
            ["short_introns_predicted_percentage","short_introns_predicted","short_introns_total"],
            ["short_introns_unpredicted_percentage","short_introns_unpredicted","short_introns_total"],
            ["long_introns_predicted_percentage","long_introns_predicted","long_introns_total"],
            ["long_introns_unpredicted_percentage","long_introns_unpredicted","long_introns_total"],
            ["start_codons_predicted_percentage","start_codons_predicted","start_codons_total"],
            ["start_codons_unpredicted_percentage","start_codons_unpredicted","start_codons_total"],
            ["stop_codons_predicted_percentage","stop_codons_predicted","stop_codons_total"],
            ["stop_codons_unpredicted_percentage","stop_codons_unpredicted","stop_codons_total"],
            ["first_exons_predicted_percentage","first_exons_predicted","first_exons_total"],
            ["first_exons_unpredicted_percentage","first_exons_unpredicted","first_exons_total"],
            ["last_exons_predicted_percentage","last_exons_predicted","last_exons_total"],
            ["last_exons_unpredicted_percentage","last_exons_unpredicted","last_exons_total"],
            ["donors_predicted_percentage","donors_predicted","donors_total"],
            ["donors_unpredicted_percentage","donors_unpredicted","donors_total"],
            ["acceptors_predicted_percentage","acceptors_predicted","acceptors_total"],
            ["acceptors_unpredicted_percentage","acceptors_unpredicted","acceptors_total"]
        ],
        100
    )

    getDivisions(
        ungroupedDf,
        generalStatistic,
        [
            ["exons_total_size_avg","exons_total_size_sum","exons_total"],
            ["exons_predicted_size_avg","exons_predicted_size_sum","exons_predicted"],
            ["exons_unpredicted_size_avg","exons_unpredicted_size_sum","exons_unpredicted"],
            ["short_exons_total_size_avg","short_exons_total_size_sum","short_exons_total"],
            ["short_exons_predicted_size_avg","short_exons_predicted_size_sum","short_exons_predicted"],
            ["short_exons_unpredicted_size_avg","short_exons_unpredicted_size_sum","short_exons_unpredicted"],
            ["long_exons_total_size_avg","long_exons_total_size_sum","long_exons_total"],
            ["long_exons_predicted_size_avg","long_exons_predicted_size_sum","long_exons_predicted"],
            ["long_exons_unpredicted_size_avg","long_exons_unpredicted_size_sum","long_exons_unpredicted"],
            ["introns_total_size_avg","introns_total_size_sum","introns_total"],
            ["introns_predicted_size_avg","introns_predicted_size_sum","introns_predicted"],
            ["introns_unpredicted_size_avg","introns_unpredicted_size_sum","introns_unpredicted"],
            ["short_introns_total_size_avg","short_introns_total_size_sum","short_introns_total"],
            ["short_introns_predicted_size_avg","short_introns_predicted_size_sum","short_introns_predicted"],
            ["short_introns_unpredicted_size_avg","short_introns_unpredicted_size_sum","short_introns_unpredicted"],
            ["long_introns_total_size_avg","long_introns_total_size_sum","long_introns_total"],
            ["long_introns_predicted_size_avg","long_introns_predicted_size_sum","long_introns_predicted"],
            ["long_introns_unpredicted_size_avg","long_introns_unpredicted_size_sum","long_introns_unpredicted"],
            ["first_exons_total_size_avg","first_exons_total_size_sum","first_exons_total"],
            ["first_exons_predicted_size_avg","first_exons_predicted_size_sum","first_exons_predicted"],
            ["first_exons_unpredicted_size_avg","first_exons_unpredicted_size_sum","first_exons_unpredicted"],
            ["last_exons_total_size_avg","last_exons_total_size_sum","last_exons_total"],
            ["last_exons_predicted_size_avg","last_exons_predicted_size_sum","last_exons_predicted"],
            ["last_exons_unpredicted_size_avg","last_exons_unpredicted_size_sum","last_exons_unpredicted"]
        ]
    )

    generalStatistic["groups"] = {}
    generalStatistic["groups"]["group_size"] = regionSize
    generalStatistic["groups"]["min_group"] = 0
    generalStatistic["groups"]["max_group"] = int(maxPos)

    for i in range(maxPos + 1):
        initializeStatisticGroup(generalStatistic,i,regionSize)
        groupDfFiltered = groupedDf.loc[groupedDf["pos"] == i, ["identifier", "value"]].copy()


        exonsQttyGroup, exonsPredictedQttyGroup, exonsUnpredictedQttyGroup, firstExonsQttyGroup, firstExonsPredictedQttyGroup, firstExonsUnpredictedQttyGroup, lastExonsQttyGroup, lastExonsPredictedQttyGroup, lastExonsUnpredictedQttyGroup, intronsQttyGroup, intronsPredictedQttyGroup, intronsUnpredictedQttyGroup = getTotals(groupDfFiltered,["exons_total","exons_predicted","exons_unpredicted","first_exons_total","first_exons_predicted","first_exons_unpredicted","last_exons_total","last_exons_predicted","last_exons_unpredicted","introns_total","introns_predicted","introns_unpredicted"])


        getIntDivisions(
            generalStatistic["groups"][str(i)]["data"],
            [
                ["exons_by_total_percentage",exonsQttyGroup,exonsQtty],
                ["exons_predicted_by_total_percentage",exonsPredictedQttyGroup,exonsPredictedQtty],
                ["exons_unpredicted_by_total_percentage",exonsUnpredictedQttyGroup,exonsUnpredictedQtty],
                ["first_exons_by_total_percentage",firstExonsQttyGroup,firstExonsQtty],
                ["first_exons_predicted_by_total_percentage",firstExonsPredictedQttyGroup,firstExonsPredictedQtty],
                ["first_exons_unpredicted_by_total_percentage",firstExonsUnpredictedQttyGroup,firstExonsUnpredictedQtty],
                ["last_exons_by_total_percentage",lastExonsQttyGroup,lastExonsQtty],
                ["last_exons_predicted_by_total_percentage",lastExonsPredictedQttyGroup,lastExonsPredictedQtty],
                ["last_exons_unpredicted_by_total_percentage",lastExonsUnpredictedQttyGroup,lastExonsUnpredictedQtty],
                ["introns_by_total_percentage",intronsQttyGroup,intronsQtty],
                ["introns_predicted_by_total_percentage",intronsPredictedQttyGroup,intronsPredictedQtty],
                ["introns_unpredicted_by_total_percentage",intronsUnpredictedQttyGroup,intronsUnpredictedQtty]
            ],
            100
        )        

        getDivisions(
            groupDfFiltered,
            generalStatistic["groups"][str(i)]["data"],
            [
                ["exons_predicted_percentage","exons_predicted","exons_total"],
                ["exons_unpredicted_percentage","exons_unpredicted","exons_total"],
                ["introns_predicted_percentage","introns_predicted","introns_total"],
                ["introns_unpredicted_percentage","introns_unpredicted","introns_total"],
                ["first_exons_predicted_percentage","first_exons_predicted","first_exons_total"],
                ["first_exons_unpredicted_percentage","first_exons_unpredicted","first_exons_total"],
                ["last_exons_predicted_percentage","last_exons_predicted","last_exons_total"],
                ["last_exons_unpredicted_percentage","last_exons_unpredicted","last_exons_total"]
            ],
            100
        )

        getDivisions(
            groupDfFiltered,
            generalStatistic["groups"][str(i)]["data"],
            [
                ["exons_total_size_avg","exons_total_size_sum","exons_total"],
                ["exons_predicted_size_avg","exons_predicted_size_sum","exons_predicted"],
                ["exons_unpredicted_size_avg","exons_unpredicted_size_sum","exons_unpredicted"],
                ["introns_total_size_avg","introns_total_size_sum","introns_total"],
                ["introns_predicted_size_avg","introns_predicted_size_sum","introns_predicted"],
                ["introns_unpredicted_size_avg","introns_unpredicted_size_sum","introns_unpredicted"],
                ["first_exons_total_size_avg","first_exons_total_size_sum","first_exons_total"],
                ["first_exons_predicted_size_avg","first_exons_predicted_size_sum","first_exons_predicted"],
                ["first_exons_unpredicted_size_avg","first_exons_unpredicted_size_sum","first_exons_unpredicted"],
                ["last_exons_total_size_avg","last_exons_total_size_sum","last_exons_total"],
                ["last_exons_predicted_size_avg","last_exons_predicted_size_sum","last_exons_predicted"],
                ["last_exons_unpredicted_size_avg","last_exons_unpredicted_size_sum","last_exons_unpredicted"]
            ]
        )

    with open(statPath,"w") as outputFile:
        json.dump(generalStatistic, outputFile, indent=4)

def statisticalAnalysis(saveFilesBasePath,extraArgs):

    executeAnalysis = not checkParam(extraArgs,"--no-statistical-analysis")[0]
    regionSizeValid, regionSizeValue = checkParam(extraArgs,"--region-size")
    regionSize = int(regionSizeValue) if regionSizeValid else 100
    writeHeaders(saveFilesBasePath,["chromosomeCSVs/recallStatistics.csv","chromosomeCSVs/precisionStatistics.csv","chromosomeCSVs/gene_transcript_predicted.csv"],["identifier,pos,value","identifier,pos,value","chromosome_identifier,gene_id,transcript_id,start_gene,end_gene,start_transcript,end_transcript,exon_qtty,gene_string,predicted,is_baseline,gene_predicted"])
    
    if executeAnalysis:
        chromosomeFolders = getChromosomeFolders(saveFilesBasePath)
        for cf in chromosomeFolders:
            generateStatisticsPerFolder(saveFilesBasePath,cf,regionSize)
            moveGeneTranscript(saveFilesBasePath,cf)
    generateGeneralStatistics(saveFilesBasePath,f"{saveFilesBasePath}chromosomeCSVs/recallStatistics.csv",f"{saveFilesBasePath}finalJsons/recallStatistics.json",regionSize,True)
    generateGeneralStatistics(saveFilesBasePath,f"{saveFilesBasePath}chromosomeCSVs/precisionStatistics.csv",f"{saveFilesBasePath}finalJsons/precisionStatistics.json",regionSize,False)
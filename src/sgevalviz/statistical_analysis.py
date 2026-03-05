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
    predictionData = [
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
    
    predictionDataDf = pd.DataFrame(predictionData, columns=["identifier", "value"])

    return predictionDataDf

def emptyDf():
    predictionData = [(pair,0) for pair in [
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

    predictionDataDf = pd.DataFrame(predictionData, columns=["identifier", "value"])

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
    baselineNucleotidesCsv = f"{chromosomeFolder}/baselineNucleotides.csv"
    recallCsv = f"{saveFilesBasePath}chromosomeCSVs/recallStatistics.csv"

    precisionChromosomeCsv = f"{chromosomeFolder}/precisionStatistics.csv"
    candidateNucleotidesCsv = f"{chromosomeFolder}/candidateNucleotides.csv"
    precisionCsv = f"{saveFilesBasePath}chromosomeCSVs/precisionStatistics.csv"

    
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


    recallData = [
        ("nucleotides_total", baselineNucleotidesSize),
        ("nucleotides_predicted", intersectionNucleotidesSize),
        ("nucleotides_unpredicted", (baselineNucleotidesSize - intersectionNucleotidesSize)),
        ("gene_total_size", recallGenesTotal),
        ("gene_ignored", recallGenesIgnored),
        ("gene_not_ignored", recallGenesNotIgnored)
    ]

    precisionData = [
        ("nucleotides_total", candidateNucleotidesSize),
        ("nucleotides_predicted", intersectionNucleotidesSize),
        ("nucleotides_unpredicted", (candidateNucleotidesSize - intersectionNucleotidesSize)),
        ("gene_total_size", precisionGenesTotal),
        ("gene_ignored", precisionGenesIgnored),
        ("gene_not_ignored", precisionGenesNotIgnored)
    ]

    recallDataDf = pd.DataFrame(recallData, columns=["identifier", "value"])
    precisionDataDf = pd.DataFrame(precisionData, columns=["identifier", "value"])


    recallDataDf.to_csv(recallChromosomeCsv, encoding='utf-8', mode="a", index=False, header=False)
    recallDataDf.to_csv(recallCsv, encoding='utf-8', mode="a", index=False, header=False)

    precisionDataDf.to_csv(precisionChromosomeCsv, encoding='utf-8', mode="a", index=False, header=False)
    precisionDataDf.to_csv(precisionCsv, encoding='utf-8', mode="a", index=False, header=False)

def generateStatisticsPerFolder(saveFilesBasePath,sf,regionSize):
    hasBaseline = os.path.exists(f"{sf}/processedBaselineFile.csv")
    hasCandidate = os.path.exists(f"{sf}/processedCandidateFile.csv")

    baselineDf = pd.read_csv(f"{sf}/processedBaselineFile.csv") if hasBaseline else None
    candidateDf = pd.read_csv(f"{sf}/processedCandidateFile.csv") if hasCandidate else None

    setNumericalStatisticsSingleFile(baselineDf,candidateDf,f"{sf}/recallStatistics.csv",f"{saveFilesBasePath}chromosomeCSVs/recallStatistics.csv",regionSize,f"{sf}/baselineNucleotides.csv")
    setNumericalStatisticsSingleFile(candidateDf,baselineDf,f"{sf}/precisionStatistics.csv",f"{saveFilesBasePath}chromosomeCSVs/precisionStatistics.csv",regionSize,f"{sf}/candidateNucleotides.csv")    
    
    writeNucleotides(saveFilesBasePath,sf)
    moveGeneTranscript(saveFilesBasePath,sf)

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
    baseDf = baseDf.groupby(["identifier"], as_index=False)["value"].sum()
    newDf = baseDf.copy()
    return newDf

def getTotals(df,cols):
    return [df.loc[df["identifier"] == col,"value"].iloc[0] for col in cols]

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

    df = breakDf(filePath)

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
        df,
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
        df,
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

def statisticalAnalysis(saveFilesBasePath,extraArgs):
    regionSizeValid, regionSizeValue = checkParam(extraArgs,"--region-size")
    regionSize = int(regionSizeValue) if regionSizeValid else 100
    writeHeaders(saveFilesBasePath,["chromosomeCSVs/recallStatistics.csv","chromosomeCSVs/precisionStatistics.csv","chromosomeCSVs/gene_transcript_predicted.csv"],["identifier,value","identifier,value","chromosome_identifier,gene_id,transcript_id,start_gene,end_gene,start_transcript,end_transcript,exon_qtty,gene_string,predicted,is_baseline,gene_predicted"])
    
    chromosomeFolders = getChromosomeFolders(saveFilesBasePath)
    for cf in chromosomeFolders:
        generateStatisticsPerFolder(saveFilesBasePath,cf,regionSize)
    generateGeneralStatistics(saveFilesBasePath,f"{saveFilesBasePath}chromosomeCSVs/recallStatistics.csv",f"{saveFilesBasePath}finalJsons/recallStatistics.json",regionSize,True)
    generateGeneralStatistics(saveFilesBasePath,f"{saveFilesBasePath}chromosomeCSVs/precisionStatistics.csv",f"{saveFilesBasePath}finalJsons/precisionStatistics.json",regionSize,False)
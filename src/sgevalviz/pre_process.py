import os
from sgevalviz.utils import *
import pandas as pd
import shutil

def trimAttribute(attributeVal):
    return attributeVal.strip("\n").strip(";").strip('"')

def getGtfExtraAttributes(extraAttributes, attributesName):
    extraAttributesList = extraAttributes.strip().split(" ")

    if len(extraAttributesList) == 1:
        return False, extraAttributesList[0]
    extraAttributesDict = {extraAttributesList[i]: trimAttribute(extraAttributesList[i+1]) for i in range(0,len(extraAttributesList),2) if extraAttributesList[i] in attributesName}
    return True, extraAttributesDict

def getLineParams(line,config,isStandardConfig):
    falseResult = [False for i in range(10)]

    if line.count("\t") != 8:
        return falseResult

    seqname, source, featureType, startPos, endPos, score, strand, frame, extraAttributes = line.split('\t')

    if featureType not in ["CDS","start_codon","stop_codon","gene","transcript"]:
        return falseResult
    
    if featureType == "gene":
        isDict , geneValue = getGtfExtraAttributes(extraAttributes, ["gene_id"])
        if isDict:
            extraAttributesDict = {"gene_id": geneValue["gene_id"], "transcript_id": None}
        else:
            extraAttributesDict = {"gene_id": geneValue, "transcript_id": None}
    elif featureType == "transcript":
        isDict , transcriptValue = getGtfExtraAttributes(extraAttributes, ["gene_id", "transcript_id"])
        if isDict:
            extraAttributesDict = {"gene_id": transcriptValue["gene_id"], "transcript_id": transcriptValue["transcript_id"]}
        else:
            extraAttributesDict = {"gene_id": transcriptValue.split(".")[0], "transcript_id": transcriptValue}
    else:
        _ , extraAttributesDict = getGtfExtraAttributes(extraAttributes, ["gene_id", "transcript_id"])

    geneId = extraAttributesDict["gene_id"]
    transcriptId = extraAttributesDict["transcript_id"]
    if config != "":
        validLine, seqname, source, featureType, startPos, endPos, score, strand, frame, geneId, transcriptId = updateLineParamsToConfig(config,isStandardConfig, seqname, source, featureType, startPos, endPos, score, strand, frame, geneId, transcriptId)

        if not validLine:
            return falseResult

    return seqname, source, featureType, startPos, endPos, score, strand, frame, geneId, transcriptId


def createRegionLine(chromosome_identifier, gene_id, transcript_id, is_exon , is_intron , is_start_codon, is_stop_codon, strand, region_start, region_end):
    is_foward_strand = ("True" if strand == "+" else "False")
    newLine = f"{chromosome_identifier},{gene_id},{transcript_id},{is_exon},{is_intron},{is_start_codon},{is_stop_codon},nan,nan,{is_foward_strand},{region_start},{region_end},nan\n"
    return newLine

def getChromosomeIdentifier(seqname,strand,splitByChromosome):
    strandName = f"{'foward' if strand == '+' else 'reverse'}_strand"
    chromosomeName = f"{seqname}__{strandName}"
    chromosomeIdentifier = (chromosomeName if splitByChromosome else strandName)

    return chromosomeIdentifier

def writeOnChromosomeFolder(input,output,chromosomeIdentifier,header):
    with open(input, 'r') as f_in, open(output, 'w') as f_out:
        f_out.write(f"{header}\n")
        for line in f_in:
            lineIdentifier, _ = line.split(",", 1)
            if lineIdentifier == chromosomeIdentifier:
                f_out.write(line)


def writeChromosomeFile(chromosomeIdentifier,baseDir,inputPath,inputPathTranscript,outputPath,transcriptGenePath):
    chromosomeNewDir = f"{baseDir}{chromosomeIdentifier}/"
    fullOutputPath = f"{chromosomeNewDir}{outputPath}"
    fullTranscriptPath = f"{chromosomeNewDir}{transcriptGenePath}" 
    os.makedirs(chromosomeNewDir,exist_ok=True)

    writeOnChromosomeFolder(inputPath,fullOutputPath,chromosomeIdentifier,"chromosome_identifier,gene_id,transcript_id,is_exon,is_intron,is_start_codon,is_stop_codon,is_first_exon,is_last_exon,is_foward_strand,region_start,region_end,predicted")
    writeOnChromosomeFolder(inputPathTranscript,fullTranscriptPath,chromosomeIdentifier,"chromosome_identifier,gene_id,transcript_id,start_gene,end_gene,start_transcript,end_transcript")

    return

def isInvalidLine(line):
    strippedLine = line.strip()
    return strippedLine == "" or strippedLine.startswith("#")

def writeSinglePreProcess(inputPath,outputPath,transcriptGenePath,splitByChromosome,extraArgs,config,isStandardConfig):
    chromosomes = set()
    geneDf = pd.DataFrame(columns=['chromosome_identifier','gene_id','start_gene','end_gene'])
    transcriptDf = pd.DataFrame(columns=['chromosome_identifier','gene_id','transcript_id','start_transcript','end_transcript'])
    with open(inputPath, 'r') as f_in, open(outputPath, 'w') as f_out:
        for line in f_in:
            if isInvalidLine(line):
                continue

            seqname, source, featureType, startPos, endPos, score, strand, frame, geneId, transcriptId = getLineParams(line,config,isStandardConfig)         

            if featureType == False:
                continue
            chromosomeIdentifier = getChromosomeIdentifier(seqname,strand,splitByChromosome)
            chromosomes.add(chromosomeIdentifier)

            if featureType == "CDS":
                exonLine = createRegionLine(chromosomeIdentifier,geneId,transcriptId,"True","False","False","False",strand,startPos,endPos)
                intronLine = createRegionLine(chromosomeIdentifier,geneId,transcriptId,"False","True","False","False",strand,str(int(endPos)+ 1),"nan")
                f_out.write(exonLine)
                f_out.write(intronLine)
            elif featureType in ['start_codon','stop_codon']:
                isStartCodon = ("True" if featureType == "start_codon" else "False")
                isStopCodon = ("True" if featureType == "stop_codon" else "False")
                startStopCodonLine = createRegionLine(chromosomeIdentifier,geneId,transcriptId,"False","False",isStartCodon,isStopCodon,strand,startPos,endPos)
                f_out.write(startStopCodonLine)
            elif featureType == "gene":
                newRow = {'chromosome_identifier':chromosomeIdentifier,'gene_id':geneId,'start_gene':startPos,'end_gene':endPos}
                geneDf = pd.concat([geneDf, pd.DataFrame([newRow])], ignore_index=True)
            elif featureType == "transcript":
                newRow = {'chromosome_identifier':chromosomeIdentifier,'gene_id':geneId,'transcript_id':transcriptId,'start_transcript':startPos,'end_transcript':endPos}
                transcriptDf = pd.concat([transcriptDf, pd.DataFrame([newRow])], ignore_index=True)

    geneTranscriptDf = pd.merge(transcriptDf,geneDf,on=['chromosome_identifier','gene_id'],how='left')
    geneTranscriptDf = geneTranscriptDf[["chromosome_identifier","gene_id","transcript_id","start_gene","end_gene","start_transcript","end_transcript"]]
    geneTranscriptDf.to_csv(transcriptGenePath, encoding='utf-8', index=False)
    return chromosomes

def preProcessFile(baseDir, inputPath,outputPath, transcriptGenePath, splitByChromosome,extraArgs,config,isStandardConfig):
    singleFileOutput = f"{baseDir}{outputPath}"
    transcriptGeneFileOutput = f"{baseDir}{transcriptGenePath}"
    chromosomeIdentifiers = writeSinglePreProcess(inputPath,singleFileOutput,transcriptGeneFileOutput,splitByChromosome,extraArgs,config,isStandardConfig)

    for chromosome in chromosomeIdentifiers:
        writeChromosomeFile(chromosome,baseDir,singleFileOutput,transcriptGeneFileOutput,outputPath,transcriptGenePath)

# Revised
def fileDefinition(saveFilesBasePath):
    basePath = f"{saveFilesBasePath}chromosomeCSVs/"

    if os.path.exists(basePath):
        shutil.rmtree(basePath)

    os.makedirs(basePath,exist_ok=True)

    finalPath = f"{saveFilesBasePath}finalJsons/"

    if os.path.exists(finalPath):
        shutil.rmtree(finalPath)

    os.makedirs(finalPath,exist_ok=True)

    processedCandidateFile = "processedCandidateFile.csv"
    processedBaselineFile = "processedBaselineFile.csv"

    transcriptAndGeneCandidateFile = "transcriptAndGeneCandidateFile.csv"
    transcriptAndGeneBaselineFile = "transcriptAndGeneBaselineFile.csv"


    return basePath, processedCandidateFile, processedBaselineFile, transcriptAndGeneCandidateFile, transcriptAndGeneBaselineFile


def preProcess(saveFilesBasePath,candidateFilePath, baselineFilePath, extraArgs):
    splitByChromosome = not checkParam(extraArgs,"--no-split")[0]
    
    standardCandidateConfig =  checkParam(extraArgs,"--candidate-config")[1] if checkParam(extraArgs,"--candidate-config")[0] else ""
    customCandidateConfig =  checkParam(extraArgs,"--custom-candidate-config")[1] if checkParam(extraArgs,"--custom-candidate-config")[0] else ""
    candidateConfig, isStandardCandidateConfig = getConfigType(standardCandidateConfig, customCandidateConfig)
    
    standardBaselineConfig =  checkParam(extraArgs,"--baseline-config")[1] if checkParam(extraArgs,"--baseline-config")[0] else ""
    customBaselineConfig =  checkParam(extraArgs,"--custom-baseline-config")[1] if checkParam(extraArgs,"--custom-baseline-config")[0] else ""
    baselineConfig, isStandardBaselineConfig = getConfigType(standardBaselineConfig, customBaselineConfig)

    basePath, processedCandidateFile, processedBaselineFile, transcriptAndGeneCandidateFile, transcriptAndGeneBaselineFile = fileDefinition(saveFilesBasePath)

    preProcessFile(basePath, candidateFilePath, processedCandidateFile, transcriptAndGeneCandidateFile, splitByChromosome,extraArgs,candidateConfig,isStandardCandidateConfig)
    preProcessFile(basePath, baselineFilePath, processedBaselineFile, transcriptAndGeneBaselineFile, splitByChromosome,extraArgs,baselineConfig,isStandardBaselineConfig)

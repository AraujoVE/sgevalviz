import pandas as pd
from sgevalviz.reader import Reader
from sgevalviz.gtf_processor import GtfProcessor

def writeProcessedLine(file, chromosomeId, gtfParams, lineType, startPos, endPos, header):
    geneId = gtfParams["gene_id"]
    transcriptId = gtfParams["transcript_id"]

    isExon = str(lineType == "exon")
    isIntron = str(lineType == "intron")
    isStartCodon = str(lineType == "start_codon")
    isStopCodon = str(lineType == "stop_codon")

    isFirstExon = "False"
    isLastExon = "False"
    isSingleExon = "False"
    isIntronRetentionExon = "False"

    isForwardStrand = str(gtfParams["strand"] == "+")

    regionStart = startPos
    regionEnd = endPos

    predicted = "False"
    genePredicted = "False"

    nucleotideSize = ""
    nucleotideList = ""

    newList = [
        chromosomeId,
        geneId,
        transcriptId,
        isExon,
        isIntron,
        isStartCodon,
        isStopCodon,
        isFirstExon,
        isLastExon,
        isSingleExon,
        isIntronRetentionExon,
        isForwardStrand,
        regionStart,
        regionEnd,
        nucleotideSize,
        nucleotideList,
        predicted,
        genePredicted
    ]

    newLine = header + ",".join(newList) + "\n"

    file.write(newLine)

def updateGeneOrTranscriptDf(df, chromosomeId, gtfParams, isGene):
    if isGene:
        newRow = {
            "chromosome_identifier": chromosomeId,
            "is_forward_strand": gtfParams['strand'] == '+',
            "gene_id": gtfParams['gene_id'],
            "start_gene": gtfParams['startPos'],
            "end_gene": gtfParams['endPos']
        }
    else:
        newRow = {
            "chromosome_identifier": chromosomeId,
            "is_forward_strand": gtfParams['strand'] == '+',
            "gene_id": gtfParams['gene_id'],
            "transcript_id": gtfParams['transcript_id'],
            "start_transcript": gtfParams['startPos'],
            "end_transcript": gtfParams['endPos']
        }

    newDf = pd.concat([df, pd.DataFrame([newRow])], ignore_index=True)

    return newDf

def writeSinglePreProcess(groupType: str, reader: Reader, processor: GtfProcessor):
    chromosomes = set()
    geneDf = pd.DataFrame(columns=reader.getGeneDfCols())
    transcriptDf = pd.DataFrame(columns=reader.getTranscriptDfCols())
    header = ','.join(reader.getProcessedDfCols()) + '\n'
    with open(reader.getInputPath(groupType), 'r') as f_in, open(reader.getSingleFilePath(groupType, True), 'w') as f_out:
        for line in f_in:
            if processor.isInvalidLine(line):
                continue

            gtfParams = processor.getGtfLineParams(line)

            if gtfParams is None:
                continue

            chromosomeId =  f"{gtfParams['seqname']}__{'forward' if gtfParams['strand'] == '+' else 'reverse'}_strand"
            chromosomes.add(chromosomeId)

            if gtfParams['featureType'] == "CDS":
                writeProcessedLine(f_out, chromosomeId, gtfParams, 'exon', gtfParams['startPos'], gtfParams['endPos'], header)
                header = ''
                writeProcessedLine(f_out, chromosomeId, gtfParams, 'intron', str(int(gtfParams['endPos']) + 1), 'nan', header)
            elif gtfParams['featureType'] in ['start_codon','stop_codon']:
                writeProcessedLine(f_out, chromosomeId, gtfParams, gtfParams['featureType'], gtfParams['startPos'], gtfParams['endPos'], header)
                header = ''
            elif gtfParams['featureType'] == "gene":
                geneDf = updateGeneOrTranscriptDf(geneDf, chromosomeId, gtfParams, True)
            elif gtfParams['featureType'] == "transcript":
                transcriptDf = updateGeneOrTranscriptDf(transcriptDf, chromosomeId, gtfParams, False)

    geneTranscriptDf = pd.merge(transcriptDf,geneDf,on=['chromosome_identifier','gene_id','is_forward_strand'],how='left')
    geneTranscriptDf.to_csv(reader.getSingleFilePath(groupType, False), encoding='utf-8', index=False)
    return chromosomes

def preProcessFile(groupType, reader):
    processor = GtfProcessor(reader, groupType)
    chromosomeIdentifiers = writeSinglePreProcess(groupType, reader, processor)

    for chromosomeId in chromosomeIdentifiers:
        reader.initializeChromosomeFolder(chromosomeId)
        reader.moveFromSingleToChromosome(groupType, chromosomeId, True)
        reader.moveFromSingleToChromosome(groupType, chromosomeId, False)

def preProcess(reader: Reader):
    preProcessFile("candidate", reader)
    preProcessFile("baseline", reader)

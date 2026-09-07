import os
import shutil
import pandas as pd
from typing import Literal

class Reader:
    def __init__(self, saveFilesBasePath, candidatePath, baselinePath, extraArgs):
        self.basePath = saveFilesBasePath.strip().rstrip("/")
        self.candidate = candidatePath.strip().rstrip("/")
        self.baseline = baselinePath.strip().rstrip("/")
        self.args = extraArgs
        self.gtfTransformation = self.getGtfTransformation()
        self.initializeFolders()

    ####################################################################################################
    ####################################################################################################
    ####################################################################################################

    ### Initialization
    # This block deals with the initialization of data

    def initializeFolders(self):
        newPreProcessFolders = [
            self.getChromosomesFolderPath(),
            self.getSingleTmpFolderPath(),
            self.getFinalDataFolderPath()
        ]

        for folder in newPreProcessFolders:
            if os.path.exists(folder):
                shutil.rmtree(folder)
            os.makedirs(folder,exist_ok=True)

        with open(self.getStatisticsFile(), "w") as f: f.write("")
        statisticsDf = pd.DataFrame(columns=["identifier","value"])
        statisticsDf.to_csv(self.getStatisticsFile(),index=False)

    def initializeChromosomeFolder(self, chromosomeId):
        chromosomePath = f"{self.getChromosomesFolderPath()}/{chromosomeId}"
        os.makedirs(chromosomePath, exist_ok=True)

    ####################################################################################################
    ####################################################################################################
    ####################################################################################################

    ### Get File Paths
    # This block deals with methods to get the paths and content of the files

    # Input paths

    def getInputPath(self, groupType):
        return self.candidate if groupType == "candidate" else self.baseline

    # Folder paths

    def getChromosomesFolderPath(self):
        return f"{self.basePath}/chromosomes"

    def getChromosomeFolderPath(self, chromosomeId):
        return f"{self.basePath}/chromosomes/{chromosomeId}"

    def getSingleTmpFolderPath(self):
        return f"{self.basePath}/singleTmp"

    def getFinalDataFolderPath(self):
        return f"{self.basePath}/finalData"

    # SingleTmp file paths

    def getSingleFilePath(self, groupType, isGtf):
        fileType = "processedGtf" if isGtf else "geneTranscript"
        return f"{self.getSingleTmpFolderPath()}/{fileType}_{groupType}.csv"

    # Chromosomes file paths

    def getChromosomePath(self, groupType, chromosomeId, isGtf):
        fileType = "processedGtf" if isGtf else "geneTranscript"
        return f"{self.getChromosomeFolderPath(chromosomeId)}/{fileType}_{groupType}.csv"

    def getDefinedChromosomePath(self, chromosomePath, groupType, fileType):
        return f"{chromosomePath}/{fileType}_{groupType}.csv"

    def getDefinedChromosomeSingleGeneStringPath(self, chromosomePath):
        return f"{chromosomePath}/singleGeneString.csv"

    def getStatisticsFile(self):
        return f"{self.basePath}/statistics.csv"

    ####################################################################################################
    ####################################################################################################
    ####################################################################################################

    ### Df columns
    # This block returns the list of columns of each Df

    def getGeneDfCols(self):
        return ['chromosome_identifier','is_forward_strand','gene_id','start_gene','end_gene']

    def getTranscriptDfCols(self):
        return ['chromosome_identifier','is_forward_strand','gene_id','transcript_id','start_transcript','end_transcript']

    def getGeneTranscriptDfCols(self):
        return ['chromosome_identifier','is_forward_strand','gene_id','transcript_id','start_gene','end_gene','start_transcript','end_transcript']

    def getProcessedDfCols(self):
        return ['chromosome_identifier', 'gene_id', 'transcript_id', 'is_exon', 'is_intron', 'is_start_codon', 'is_stop_codon', 'is_first_exon', 'is_last_exon', 'is_single_exon', 'is_intron_retention_exon', 'is_forward_strand', 'region_start', 'region_end', 'nucleotide_size', 'nucleotide_list','predicted', 'gene_predicted']

    def getGeneStringDfCols(self):
        return ["chromosome_identifier","gene_id", "transcript_id","start_gene","end_gene","start_transcript","end_transcript","min_pos","max_pos","strand","exon_qtty","intron_retention_qtty","gene_string","is_forward_strand","is_baseline","same_strand","predicted","gene_predicted"]

    def getStatisticsDf(self):
        return ["identifier","value"]        

    ####################################################################################################
    ####################################################################################################
    ####################################################################################################

    ### Move files
    # This block moves data from one path to another

    def moveFromSingleToChromosome(self, groupType, chromosomeId, isGtf):
        input = self.getSingleFilePath(groupType, isGtf)
        output = self.getChromosomePath(groupType, chromosomeId, isGtf)

        with open(input, 'r') as f_in, open(output, 'w') as f_out:
            for i, line in enumerate(f_in):
                if i == 0:
                    f_out.write(line)
                else:
                    lineIdentifier, _ = line.split(",", 1)
                    if lineIdentifier == chromosomeId:
                        f_out.write(line)

    ####################################################################################################
    ####################################################################################################
    ####################################################################################################

    ### GTF Transformation
    # This block deals with methods used when the gtf data should pass for a transformation that uses a config file. The config can be of a standard one or a custom file.

    def getGtfTransformation(self):
        standardCandidateConfig =  self.getParam("candidate-config")
        customCandidateConfig =  self.getParam("custom-candidate-config")
        standardBaselineConfig =  self.getParam("baseline-config")
        customBaselineConfig =  self.getParam("custom-baseline-config")

        gtfTransformation = {
            "candidate": {
                "name": customCandidateConfig or standardCandidateConfig,
                "type": "custom" if customCandidateConfig else "standard",
                "transformation": (customCandidateConfig or standardCandidateConfig) != ""
            },
            "baseline": {
                "name": customBaselineConfig or standardBaselineConfig,
                "type": "custom" if customBaselineConfig else "standard",
                "transformation": (customBaselineConfig or standardBaselineConfig) != ""
            }
        }

        return gtfTransformation

    def getGtfTransformationName(self, groupType):
        return self.gtfTransformation[groupType]["name"]

    def getGtfTransformationType(self, groupType):
        return self.gtfTransformation[groupType]["type"]


    ####################################################################################################
    ####################################################################################################
    ####################################################################################################

    ### Argument verification
    # This block deals with methods used to check the extra parameters passed to the function

    def getParam(self, param, defaultValue=""):
        for arg in self.args:
            key, sep, value = arg.partition("=")
            if key == f"--{param}":
                return value if sep else defaultValue
        return defaultValue

    ####################################################################################################
    ####################################################################################################
    ####################################################################################################

    ### Others
    # Others

    def getChromosomeFoldersList(self):
        folder = self.getChromosomesFolderPath()
        subfolders = [f.path.rstrip("/") for f in os.scandir(folder) if f.is_dir()]
        
        return subfolders

    def hasGtfFile(self, chromosomePath, groupType):
        filePath = self.getDefinedChromosomePath(chromosomePath, groupType, "processedGtf")

        return os.path.isfile(filePath)

    ####################################################################################################
    ####################################################################################################
    ####################################################################################################

    ### Statistics Writing
    # Statistics Writing

    def setFinalResults(self):
        statisticsPath = self.getStatisticsFile()
        statisticsDf = pd.read_csv(statisticsPath)
        groupedDf = statisticsDf.groupby("identifier")["value"].sum().reset_index(drop=True)

        groupedDf.to_csv(statisticsPath, index=False)


    def updateReference(self, baseName, partialValueBase, totalValue, multiplier = 1):
        statisticsPath = self.getStatisticsFile()
        statisticsDf = pd.read_csv(statisticsPath)
        partialValue = partialValueBase * multiplier
        newRows = pd.DataFrame([{"identifier": f"partial___{baseName}", "value": partialValue}, {"identifier": f"total___{baseName}", "value": totalValue}])
        statisticsDf = pd.concat([statisticsDf, newRows], ignore_index=True)
        statisticsDf.to_csv(statisticsPath, index=False)


    def statisticsUpdate__Strand__Value(
            self,
            firstLevel: Literal["forward", "general", "reverse"],
            lastLevel: Literal["reference_gene_unpredicted-percentage", "no_prediction_for_reference_on_chromosome_strand-percentage"],
            partialValue: int,
            totalValue:  int
        ):
        baseName = f"{firstLevel}__{lastLevel}"
        multiplier = 100 if lastLevel.endswith("-percentage") else 1
        self.updateReference(baseName, partialValue, totalValue, multiplier)
        return

    def statisticsUpdate__StrandRefGenePartPred__Value(
            self,
            firstLevel: Literal["forward", "general", "reverse"],
            lastLevel: Literal["selected_model_transcript_on_same_frame_as_selected_reference_transcript-percentage", "ratio_of_number_of_exons_in_model_selected_transcript_per_reference_selected_transcript-average", "ratio_of_number_of_nucleotides_in_model_selected_transcript_per_reference_selected_transcript-average"],
            partialValue: int,
            totalValue:  int
        ):
        baseName = f"{firstLevel}__reference_gene_partially_predicted__{lastLevel}"
        multiplier = 100 if lastLevel.endswith("-percentage") else 1
        self.updateReference(baseName, partialValue, totalValue, multiplier)
        return

    def statisticsUpdate__StrandRefGenePredHasIntronRetExonInModel__Value(
            self,
            firstLevel: Literal["forward", "general", "reverse"],
            lastLevel: Literal["reference_gene_partially_predicted", "reference_gene_predicted", "reference_gene_totally_predicted"],
            partialValue: int,
            totalValue:  int
        ):
        baseName = f"{firstLevel}__{lastLevel}__has_intron_retention_exon_in_reference__selected_reference_transcript_is_the_transcript_with_most_intron_retention_exons-percentage"
        multiplier = 100 if lastLevel.endswith("-percentage") else 1
        self.updateReference(baseName, partialValue, totalValue, multiplier)
        return

    def statisticsUpdate__Strand_RefGenePred_ModelTransc__Value(
            self,
            firstLevel: Literal["forward", "general", "reverse"],
            secondLevel: Literal["reference_gene_partially_predicted", "reference_gene_predicted", "reference_gene_totally_predicted"],
            thirdLevel: Literal["selected_model_transcript", "average_model_transcript"],
            lastLevel: Literal["average_size_of_exons_per_transcript-average", "number_of_exons_per_transcript-average", "average_size_of_introns_per_transcript-average"],
            partialValue: int,
            totalValue:  int
        ):
        baseName = f"{firstLevel}__{secondLevel}__{thirdLevel}__{lastLevel}"
        multiplier = 100 if lastLevel.endswith("-percentage") else 1
        self.updateReference(baseName, partialValue, totalValue, multiplier)
        return

    def statisticsUpdate__StrandRefGenePredRecallPrecision__Value(
            self,
            firstLevel: Literal["forward", "general", "reverse"],
            lastLevel: Literal["gene_predicted_recall", "gene_predicted_precision"],
            partialValue: int,
            totalValue:  int
        ):
        baseName = f"{firstLevel}__{lastLevel}__single_exon_occurance-percentage"
        multiplier = 100 if lastLevel.endswith("-percentage") else 1
        self.updateReference(baseName, partialValue, totalValue, multiplier)
        return

    def statisticsUpdate__Strand_RefGenePredRecallPrecision_ExonNumberPerTranscript__Value(
            self,
            firstLevel: Literal["forward", "general", "reverse"],
            secondLevel: Literal["gene_predicted_recall", "gene_predicted_precision"],
            thirdLevel: Literal["single_exon_selected_model_transcript", "any_quantity_exon_selected_model_transcript", "multiple_exon_selected_model_transcript"],
            lastLevel: Literal["totally_predicted_genes_prediction-percentage", "nucleotide_prediction-percentage", "start_codon_prediction-percentage", "stop_codon_prediction-percentage", "same_time_start_codon_and_stop_codon_prediction-percentage", "average_intron_prediction-percentage", "average_exon_prediction-percentage", "first_exon_prediction-percentage", "last_exon_prediction-percentage", "same_time_first_exon_and_last_exon_prediction-percentage", "average_donnor_prediction-percentage", "average_acceptor_prediction-percentage"],
            partialValue: int,
            totalValue:  int
        ):
        baseName = f"{firstLevel}__{secondLevel}__{thirdLevel}__{lastLevel}"
        multiplier = 100 if lastLevel.endswith("-percentage") else 1
        self.updateReference(baseName, partialValue, totalValue, multiplier)
        return

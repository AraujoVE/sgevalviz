import os
import shutil
import pandas as pd

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
        newPlotFolders = [
            self.getPlotsFolderPath("candidate"),
            self.getPlotsFolderPath("baseline")
        ]

        if not self.hasParam("no-pre-process"):
            for folder in newPreProcessFolders:
                if os.path.exists(folder):
                    shutil.rmtree(folder)
                os.makedirs(folder,exist_ok=True)

        if not self.hasParam("no-plot"):
            for folder in newPreProcessFolders:
                if os.path.exists(folder):
                    shutil.rmtree(folder)
                os.makedirs(folder,exist_ok=True)

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

    def getPlotsFolderPath(self, groupType):
        return f"{self.basePath}/{"Precision" if groupType == 'candidate' else "Recall"}Plots"

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

    def getSingleGeneStringPath(self):
        return f"{self.getFinalDataFolderPath()}/singleGeneString.csv"

    def getDefinedChromosomeStatistics(self, chromosomePath, groupType):
        statisticType = "recall" if groupType == "baseline" else "precision"
        return f"{chromosomePath}/{statisticType}Statistics.csv"

    def getFinalStatistics(self, groupType):
        statisticType = "recall" if groupType == "baseline" else "precision"
        return f"{self.getFinalDataFolderPath()}/{statisticType}Statistics.csv"

    def getFinalStatisticsJson(self, groupType):
        statisticType = "recall" if groupType == "baseline" else "precision"
        return f"{self.getFinalDataFolderPath()}/{statisticType}Statistics.json"


    def getDefinedChromosomeNucleotides(self, chromosomePath, groupType):
        return f"{chromosomePath}/Nucleotides_{groupType}.csv"

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
        return ['chromosome_identifier', 'gene_id', 'transcript_id', 'is_exon', 'is_intron', 'is_start_codon', 'is_stop_codon', 'is_first_exon', 'is_last_exon', 'is_intron_retention_exon', 'is_forward_strand', 'region_start', 'region_end', 'predicted', 'gene_predicted']

    def getGeneStringDfCols(self):
        return ["chromosome_identifier","gene_id", "transcript_id","start_gene","end_gene","start_transcript","end_transcript","min_pos","max_pos","strand","exon_qtty","intron_retention_qtty","gene_string","is_forward_strand","is_baseline","same_strand","predicted","gene_predicted"]

    def getStatisticsDf(self):
        return ["identifier","value"]

    def getNucleotidesDf(self):
        return []

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
    
    def needsGtfTransformation(self, groupType):
        return self.gtfTransformation[groupType]["transformation"]

    def getGtfTransformationName(self, groupType):
        return self.gtfTransformation[groupType]["name"]

    def getGtfTransformationType(self, groupType):
        return self.gtfTransformation[groupType]["type"]


    ####################################################################################################
    ####################################################################################################
    ####################################################################################################

    ### Argument verification
    # This block deals with methods used to check the extra parameters passed to the function

    def hasParam(self, param):
        return any(arg.split("=", 1)[0] == f"--{param}" for arg in self.args)

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

    def getProcessedData(self, chromosomePath, groupType):
        filePath = self.getDefinedChromosomePath(chromosomePath, groupType, "processedGtf")
        hasFile = os.path.exists(filePath)

        return pd.read_csv(filePath) if hasFile else None
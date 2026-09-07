from sgevalviz.reader import Reader
from itertools import pairwise
import json
import importlib.resources as resources
import os
import re

class GtfProcessor:
    def __init__(self, reader: Reader, groupType: str):
        self.reader = reader
        self.groupType = groupType
        self.defaultGtfParams = ["seqname", "source", "featureType", "startPos", "endPos", "score", "strand", "frame", "extraAttributes"]
        self.validFeatureTypes = ["CDS", "start_codon", "stop_codon", "gene", "transcript"]
        self.setTransformation()


    ####################################################################################################
    ####################################################################################################
    ####################################################################################################

    ### Initialization
    # This block deals with the initialization of data

    def setTransformation(self):
        standardConfig =  self.reader.getParam(f"{self.groupType}-config")
        customConfig =  self.reader.getParam(f"custom-{self.groupType}-config")

        self.name = customConfig or standardConfig
        self.isCustomTransformation = True if customConfig else False
        self.transformation = (customConfig or standardConfig) != ""

    ####################################################################################################
    ####################################################################################################
    ####################################################################################################

    ### Get gtf line parameters
    # This block returns the internal data

    def getFirstGtfParams(self, line: str):
        splittedLine = line.split("\t")
        gtfParams = [(k, v) for k, v in zip(self.defaultGtfParams, splittedLine)]

        return dict(gtfParams)

    def extractGeneAndTranscriptId(self, gtfParams: dict[str, str]):
        extraAttributesList = " ".join(gtfParams["extraAttributes"].strip().split()).split()
        featureType = gtfParams["featureType"]
        extraAttribute = None

        if len(extraAttributesList) == 1:
            extraAttribute = extraAttributesList[0]
        else:
            values = [v.strip("\n").strip(";").strip('"') for v in extraAttributesList[1::2] if v != ""]
            extraAttributesDict = {k: v for k, v in zip(extraAttributesList[0::2], values)}


        if featureType == "gene" and extraAttribute:
            geneId, transcriptId = (extraAttribute, None)
        elif featureType == "gene":
            geneId, transcriptId = (extraAttributesDict["gene_id"], None)
        elif featureType == "transcript" and extraAttribute:
            geneId, transcriptId = (extraAttribute.split(".")[0], extraAttribute)
        else:
            geneId, transcriptId = (extraAttributesDict["gene_id"], extraAttributesDict["transcript_id"])

        return geneId, transcriptId

    def loadConfig(self):
        configName = self.reader.getGtfTransformationName(self.groupType)
        isStandardConfig = self.reader.getGtfTransformationType(self.groupType) == "standard"
        try:
            if isStandardConfig:
                json_str = resources.read_text("sgevalviz.configs", f"{configName}.json", encoding="utf-8")
            else:
                if not os.path.isabs(configName):
                    raise ValueError(f"Expected an absolute path for config, got: {configName}")
                if not os.path.isfile(configName):
                    raise FileNotFoundError(f"Config file not found: {configName}")
                with open(configName, "r", encoding="utf-8") as f:
                    json_str = f.read()
            return json.loads(json_str)
        except FileNotFoundError:
            raise FileNotFoundError(f"Config file '{configName}.json' not found in sgevalviz/configs/")
        except json.JSONDecodeError as e:
            raise ValueError(f"Invalid JSON in config file '{configName}.json': {e}")


    def updatedParam(self, paramName, paramValue, jsonData):
        if paramName not in jsonData: return paramValue
        
        paramData = jsonData[paramName]

        if ("from_pattern" in paramData) ^ ("to_pattern" in paramData):
            raise ValueError(f"Invalid config. The config for the parameter {paramName} has only one of 'from_pattern' and 'to_pattern' in the config.")

        newValue = paramValue
        if "from_pattern" in paramData:
            pattern = re.compile(paramData["from_pattern"])
            newValue = re.sub(pattern,paramData["to_pattern"],paramValue)


        if "valid_pattern" in paramData:
            pattern = re.compile(paramData["valid_pattern"])
            validated = re.fullmatch(pattern, newValue) is not None
            if not validated:
                raise ValueError(f"Parameter {paramName} didn't pass the validation")

        return newValue


    def updateGtfParams(self, gtfParams: dict[str, str]):
        if not self.transformation: return gtfParams

        jsonData = self.loadConfig()

        for k, v in gtfParams.items():
            gtfParams[k] = self.updatedParam(k, v, jsonData)

        return gtfParams


    def getGtfLineParams(self, line: str):
        if line.count("\t") != 8: return None

        gtfParams = self.getFirstGtfParams(line)

        if gtfParams["featureType"] not in self.validFeatureTypes: return None

        gtfParams["gene_id"], gtfParams["transcript_id"] = self.extractGeneAndTranscriptId(gtfParams)

        gtfParams = self.updateGtfParams(gtfParams)

        return gtfParams


    ####################################################################################################
    ####################################################################################################
    ####################################################################################################

    ### Utils
    # This block have usefull methods

    def isInvalidLine(self, line: str):
        strippedLine = line.strip()
        return strippedLine == "" or strippedLine.startswith("#")
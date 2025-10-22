import json
import importlib.resources as resources
import re

parameters = {
    "--no-pre-process": "zero",
    "--no-split": "zero",
    "--candidate-config": "one",
    "--baseline-config": "one",
    "--region-size": "one"
}

parameterType = {
    "zero": "no arguments (no '=' after the argument)",
    "one": "a single argument (just a value after the argument)",
    "multiple": "multiple values (after the '=' put multiple values separated by ';', if there's just one value, just put ';' at the end)"
}

def getArgType(arg):
    if "=" not in arg:
        return "zero"
    if ";" not in arg:
        return "one"
    return "multiple"

def validateParams(args):
    strippedParams = [arg.split("=")[0] for arg in args]
    paramsDontRepeat = len(strippedParams) == len(set(strippedParams))
    if not paramsDontRepeat:
        errorMessage = "Repeated Parameters"
        return False, errorMessage

    for arg in args:
        if not arg.startswith("--"):
            errorMessage = f"Argument {arg} does not start with '--'"
            return False, errorMessage
        
        if arg.count("=") > 1:
            errorMessage = f"Argument {arg} has more than one '='"
            return False, errorMessage

        argKey = arg.split("=")[0]

        if argKey not in parameters:
            errorMessage = f"Argument of type {argKey} is not valid"
            return False, errorMessage

        argType = getArgType(arg)
        
        if parameters[argKey] != argType:
            errorMessage = f"Argument of type {argKey} requires {parameterType[parameters[argKey]]}, but you put {argType}"
            return False, errorMessage

    return True, ""




def checkParam(args,chosenArg):
    validArgs = [arg for arg in args if arg.startswith(chosenArg)]

    if validArgs == []:
        return False, None

    validArg = validArgs[0]
    if "=" not in validArg:
        return True, None
    
    validArgValue = validArg.split("=")[1]
    if ";" not in validArgValue:
        return True, validArgValue
    
    validArgValues = [("" if a == ":NoArg:" else a) for a in validArgValue.split(";") if a != ""]
    return True, validArgValues



def updatedParam(paramVar, paramName, jsonData):
    if paramName not in jsonData:
        return True, paramVar
    
    paramData = jsonData[paramName]

    if ("from_pattern" in paramData) ^ ("to_pattern" in paramData):
        return  False, "Invalid config. The config has only one of 'from_pattern' and 'to_pattern' in th config."

    returnVar = paramVar
    if "from_pattern" in paramData:
        pattern = re.compile(paramData["from_pattern"])
        returnVar = re.sub(pattern,paramData["to_pattern"],paramVar)


    if "valid_pattern" in paramData:
        pattern = re.compile(paramData["valid_pattern"])
        validated = re.fullmatch(pattern,paramVar) is not None
        if not validated:
            return False, "Invalid param"

    return True, returnVar            


from importlib import resources
import json

def loadConfig(config):
    try:
        json_str = resources.read_text("sgevalviz.configs", f"{config}.json", encoding="utf-8")
        return json.loads(json_str)
    except FileNotFoundError:
        raise FileNotFoundError(f"Config file '{config}.json' not found in sgevalviz/configs/")
    except json.JSONDecodeError as e:
        raise ValueError(f"Invalid JSON in config file '{config}.json': {e}")



def updateLineParamsToConfig(config, seqname, source, featureType, startPos, endPos, score, strand, frame, geneId, transcriptId):
    # with open(f'configs/{config}.json', 'r') as f:
    #     jsonData = json.load(f)

    jsonData = loadConfig(config)

    falseReturn = [False for i in range(11)]
    
    continueCheck, seqname = updatedParam(seqname,"seqname", jsonData)
    if not continueCheck:
        return falseReturn
    
    continueCheck, source = updatedParam(source,"source", jsonData)
    if not continueCheck:
        return falseReturn
    
    continueCheck, featureType = updatedParam(featureType,"featureType", jsonData)
    if not continueCheck:
        return falseReturn
    
    continueCheck, startPos = updatedParam(startPos,"startPos", jsonData)
    if not continueCheck:
        return falseReturn
    
    continueCheck, endPos = updatedParam(endPos,"endPos", jsonData)
    if not continueCheck:
        return falseReturn
    
    continueCheck, score = updatedParam(score,"score", jsonData)
    if not continueCheck:
        return falseReturn
    
    continueCheck, strand = updatedParam(strand,"strand", jsonData)
    if not continueCheck:
        return falseReturn
    
    continueCheck, frame = updatedParam(frame,"frame", jsonData)
    if not continueCheck:
        return falseReturn
    
    continueCheck, geneId = updatedParam(geneId,"geneId", jsonData)
    if not continueCheck:
        return falseReturn
    
    continueCheck, transcriptId = updatedParam(transcriptId,"transcriptId", jsonData)
    if not continueCheck:
        return falseReturn

    return True, seqname, source, featureType, startPos, endPos, score, strand, frame, geneId, transcriptId
    

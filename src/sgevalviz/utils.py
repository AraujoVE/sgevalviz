import json
import importlib.resources as resources
import re
import os

parameters = {
    "--no-pre-process": "zero",
    "--no-plot": "zero",
    "--no-split": "zero",
    "--candidate-config": "one",
    "--baseline-config": "one",
    "--custom-candidate-config": "one",
    "--custom-baseline-config": "one",
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

def validateInputs(argv):    
    validStatus = False
    resultMsg = ""
    if len(argv) < 4:
        resultMsg = "Invalid number of arguments, should be: sgevalviz <savePath> <candidatePath> <baselinePath> [options]" 
        return validStatus, resultMsg, None, None, None

    saveFilesBasePath = argv[1]
    if not os.path.isdir(saveFilesBasePath):
        resultMsg = "Invalid argument: the first argument must be a folder"
        return validStatus, resultMsg, None, None, None

    saveFilesBasePath += "/" if saveFilesBasePath[-1] != "/" else ""

    candidatePath = argv[2]
    if not os.path.isfile(candidatePath):
        resultMsg = "Invalid argument: the second argument must be a file"
        return validStatus, resultMsg, None, None, None


    baselinePath = argv[3]
    if not os.path.isfile(baselinePath):
        resultMsg = "Invalid argument: the third argument must be a file"
        return validStatus, resultMsg, None, None, None

    validStatus = True
    return validStatus, resultMsg, saveFilesBasePath, candidatePath, baselinePath
import sys
from sgevalviz.pre_process import preProcess
from sgevalviz.fill_data import fillData
from sgevalviz.statistical_analysis import statisticalAnalysis
from sgevalviz.utils import validateParams, validateInputs, checkParam
from sgevalviz.plot import plot


def run(argv=None):
    if argv is None:
        argv = sys.argv

    if len(argv) < 4:
        print("Usage: sgevalviz <savePath> <candidatePath> <baselinePath> [options]")
        sys.exit(1)

    validStatus, resultMsg, saveFilesBasePath, candidatePath, baselinePath = validateInputs(argv)
    if validStatus == False:
        print(resultMsg)
        sys.exit(1)

    extraArgs = argv[4:]
    validArgs, errorMessage = validateParams(extraArgs)
    if not validArgs:
        print(errorMessage)
        sys.exit(1)

    if not checkParam(extraArgs,"--no-pre-process")[0]:
        preProcess(saveFilesBasePath,candidatePath, baselinePath, extraArgs)
        fillData(saveFilesBasePath,extraArgs)
        statisticalAnalysis(saveFilesBasePath,extraArgs)
    if not checkParam(extraArgs, "--no-plot")[0]:
        plot(saveFilesBasePath)


def main():
    run(sys.argv)

if __name__ == "__main__":
    main()
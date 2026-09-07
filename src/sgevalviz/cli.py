import sys
from sgevalviz.pre_process import preProcess
from sgevalviz.fill_data import fillData
from sgevalviz.utils import validateParams, validateInputs
from sgevalviz.reader import Reader


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

    reader = Reader(saveFilesBasePath, candidatePath, baselinePath, extraArgs)

    preProcess(reader)
    fillData(reader)


def main():
    run(sys.argv)

if __name__ == "__main__":
    main()
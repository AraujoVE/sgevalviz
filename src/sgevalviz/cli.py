import sys
from sgevalviz.pre_process import preProcess
from sgevalviz.fill_data import fillData
from sgevalviz.statistical_analysis import statisticalAnalysis
from sgevalviz.utils import validateParams
from sgevalviz.plot import plot


def run(argv=None):
    if argv is None:
        argv = sys.argv

    if len(argv) < 4:
        print("Usage: sgevalviz <savePath> <candidatePath> <baselinePath> [options]")
        sys.exit(1)

    saveFilesBasePath = argv[1]
    candidatePath = argv[2]
    baselinePath = argv[3]
    extraArgs = argv[4:]

    validArgs, errorMessage = validateParams(extraArgs)
    if not validArgs:
        print(errorMessage)
        sys.exit(1)

    preProcess(saveFilesBasePath,candidatePath, baselinePath, extraArgs)
    fillData(saveFilesBasePath,extraArgs)
    statisticalAnalysis(saveFilesBasePath,extraArgs)
    plot(saveFilesBasePath)


def main():
    run(sys.argv)

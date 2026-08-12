import matplotlib.pyplot as plt
import json
import os
import shutil
from sgevalviz.reader import Reader




def plot(reader: Reader):
    plotGroup(saveFilesBasePath,"Precision",f"{saveFilesBasePath}finalJsons/precisionStatistics.json")
    plotGroup(saveFilesBasePath,"Recall",f"{saveFilesBasePath}finalJsons/recallStatistics.json")



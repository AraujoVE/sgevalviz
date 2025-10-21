import matplotlib.pyplot as plt
import json
import os
import shutil

percentageConfig = [
    [
        ("nucleotides_predicted_percentage","Nucleotides"),
        ("genes_predicted_percentage","Genes"),
        ("start_codons_predicted_percentage","Start Codon"),
        ("stop_codons_predicted_percentage","Stop Codon"),
        ("genes_not_ignored_percentage", "Genes Not Ignored")
    ],
    [
        ("exons_predicted_percentage","General"),
        ("short_exons_predicted_percentage","<= 300 nt"),
        ("long_exons_predicted_percentage","> 300 nt"),
        ("first_exons_predicted_percentage","First"),
        ("last_exons_predicted_percentage","Last")
    ],
    [
        ("introns_predicted_percentage","General"),
        ("short_introns_predicted_percentage","<= 300 nt"),
        ("long_introns_predicted_percentage","> 300 nt"),
        ("donors_predicted_percentage","Donors"),
        ("acceptors_predicted_percentage","Acceptors")
    ]
]
percentageTitles = ("General Prediction","Exon Prediction","Intron Prediction")

groupPercentageConfig = [
    [
        ("exons_predicted_percentage","Exons"),
        ("first_exons_predicted_percentage","First Exons"),
        ("last_exons_predicted_percentage","Last Exons"),
        ("introns_predicted_percentage","Introns")
    ],
    [
        ("exons_by_total_percentage","General"),
        ("exons_predicted_by_total_percentage","Predicted"),
        ("exons_unpredicted_by_total_percentage","Unpredicted")
    ],
    [
        ("first_exons_by_total_percentage","General"),
        ("first_exons_predicted_by_total_percentage","Predicted"),
        ("first_exons_unpredicted_by_total_percentage","Unpredicted")
    ],
    [
        ("last_exons_by_total_percentage","General"),
        ("last_exons_predicted_by_total_percentage","Predicted"),
        ("last_exons_unpredicted_by_total_percentage","Unpredicted")
    ],
    [
        ("introns_by_total_percentage","General"),
        ("introns_predicted_by_total_percentage","Predicted"),
        ("introns_unpredicted_by_total_percentage","Unpredicted")
    ]
]


"exons_predicted_by_total_percentage"
groupPercentageTitles = ("Exon/Intron Prediction - Group","Exon Occurance Percentage","First Exon Occurance Percentage","Last Exon Occurance Percentage","Intron Occurance Percentage")


exonSizeConfig = [
    [
        ("exons_total_size_avg","All"),
        ("exons_predicted_size_avg","Predicted"),
        ("exons_unpredicted_size_avg","Unpredicted")
    ],
    [
        ("short_exons_total_size_avg","All"),
        ("short_exons_predicted_size_avg","Predicted"),
        ("short_exons_unpredicted_size_avg","Unpredicted")
    ],
    [
        ("long_exons_total_size_avg","All"),
        ("long_exons_predicted_size_avg","Predicted"),
        ("long_exons_unpredicted_size_avg","Unpredicted")
    ],
    [
        ("first_exons_total_size_avg","All"),
        ("first_exons_predicted_size_avg","Predicted"),
        ("first_exons_unpredicted_size_avg","Unpredicted")
    ],
    [
        ("last_exons_total_size_avg","All"),
        ("last_exons_predicted_size_avg","Predicted"),
        ("last_exons_unpredicted_size_avg","Unpredicted")
    ]
]
exonSizeTitles = ["General Exon Size","Short Exon Size","Long Exon Size","General First Exon Size","General Last Exon Size"]

groupSizeConfig = [
    [
        ("exons_total_size_avg","All"),
        ("exons_predicted_size_avg","Predicted"),
        ("exons_unpredicted_size_avg","Unpredicted")
    ],
    [
        ("first_exons_total_size_avg","All"),
        ("first_exons_predicted_size_avg","Predicted"),
        ("first_exons_unpredicted_size_avg","Unpredicted")
    ],
    [
        ("last_exons_total_size_avg","All"),
        ("last_exons_predicted_size_avg","Predicted"),
        ("last_exons_unpredicted_size_avg","Unpredicted")
    ],
    [
        ("introns_total_size_avg","All"),
        ("introns_predicted_size_avg","Predicted"),
        ("introns_unpredicted_size_avg","Unpredicted")
    ]
]
groupSizeTitles = ["General Exon Size - Group","General First Exon Size - Group","General Last Exon Size - Group","General Intron Size - Group"]

groupIntronSizeConfig = [
    [
        ("introns_total_size_avg","All"),
        ("introns_predicted_size_avg","Predicted"),
        ("introns_unpredicted_size_avg","Unpredicted")
    ],
    [
        ("short_introns_total_size_avg","All"),
        ("short_introns_predicted_size_avg","Predicted"),
        ("short_introns_unpredicted_size_avg","Unpredicted")
    ],
    [
        ("long_introns_total_size_avg","All"),
        ("long_introns_predicted_size_avg","Predicted"),
        ("long_introns_unpredicted_size_avg","Unpredicted")
    ]
]
groupIntronSizeTitles = ["General Intron - Group","Short Intron Size - Group","Long Intron Size - Group"]


intronSizeConfig = [
    [
        ("introns_total_size_avg","All"),
        ("introns_predicted_size_avg","Predicted"),
        ("introns_unpredicted_size_avg","Unpredicted")
    ],
    [
        ("short_introns_total_size_avg","All"),
        ("short_introns_predicted_size_avg","Predicted"),
        ("short_introns_unpredicted_size_avg","Unpredicted")
    ],
    [
        ("long_introns_total_size_avg","All"),
        ("long_introns_predicted_size_avg","Predicted"),
        ("long_introns_unpredicted_size_avg","Unpredicted")
    ]
]
intronSizeTitles = ["General Intron","Short Intron Size","Long Intron Size"]


avgExonsInGenesConfig = [
    [
        ("exons_in_predicted_genes_predicted_transcripts","Predicted Genes, Predicted Transcripts"),
        ("exons_in_predicted_genes_unpredicted_transcripts","predicted Genes, Unpredicted Transcripts"),
        ("exons_in_unpredicted_genes","Unpredicted Genes")
    ]
]
avgExonsInGenesTitle = ["Average Number of Exons in Genes"]

def baseDicts(qtty,keys):
    data = []
    for i in range(qtty):
        newDict = {}
        for j in keys:
            newDict[j] = []
        data.append(newDict)
    
    return data

def addTitles(data,titles):
    for i in range(len(data)):
        data[i]["title"] = titles[i]

def addBar(fullDict,dictKey,category,data,pos):
    dictValue = fullDict[dictKey]
    dictValue = -1 if dictValue == None else dictValue
    data[pos]["values"].append(dictValue)
    data[pos]["categories"].append(category)


def addLine(fullDict,dictKey,category,data,pos):
    divideByGroup = ("size_avg" in dictKey)
    groupSize = fullDict["groups"]["group_size"]
    minGroup = fullDict["groups"]["min_group"]
    maxGroup = fullDict["groups"]["max_group"]

    xValues = [i*groupSize for i in range(minGroup,maxGroup+1)]
    yValues = []

    for i in range(minGroup,maxGroup+1):
        dictValue = fullDict["groups"][str(i)]["data"][dictKey]
        if divideByGroup and (dictValue is not None):
            minValue = groupSize*(minGroup+i)
            dictValue = (dictValue - minValue)/groupSize
        yValues.append(dictValue)

    data[pos]["valuesX"].append(xValues)
    data[pos]["valuesY"].append(yValues)
    data[pos]["categories"].append(category)

def addColors(data,keyName):
    colors = [
        "#4E79A7",  # muted blue
        "#F28E2B",  # orange
        "#E15759",  # red
        "#76B7B2",  # teal
        "#59A14F",  # green
        "#EDC948",  # yellow
        "#B07AA1",  # purple
        "#FF9DA7",  # pink
        "#9C755F",  # brown
        "#BAB0AC",  # gray
        "#1F77B4",  # classic blue
        "#FF7F0E",  # vivid orange
        "#2CA02C",  # vivid green
        "#D62728",  # vivid red
        "#9467BD",  # vivid purple
        "#8C564B",  # muted brown
        "#E377C2",  # magenta
        "#7F7F7F",  # neutral gray
        "#BCBD22",  # lime
        "#17BECF"   # cyan
    ]

    for i in range(len(data)):
        valueSize = len(data[i][keyName])
        data[i]["colors"] = colors[:valueSize]

def plotBarData(data, axs, yLable):
    
    for i, d in enumerate(data):
        correctedValues = [(0 if v == -1 else v) for v in d["values"]]
        bars = axs[i].bar(d["categories"], correctedValues, color=d["colors"])
        middleHeight = (max(correctedValues) - min(correctedValues))/2

        for bar, value in zip(bars, d["values"]):
            if value == -1:
                label = "N/A"
                color = "black"
            else:
                label = f"{value:.2f}"
                color = "white"

            axs[i].text(
                bar.get_x() + bar.get_width() / 2,           # x-position (center)
                bar.get_height() / 2 if value != -1 else middleHeight,  # y-position
                label,
                ha='center',
                va='center',
                fontsize=9,
                fontweight='bold',
                color=color
            )

        axs[i].set_title(d["title"], fontsize=12, fontweight='bold')
        axs[i].set_ylabel(yLable)
        axs[i].grid(axis='y', linestyle='--', alpha=0.5)

def plotLineData(data, axs, yLable):
    for i, d in enumerate(data):
        x = d["valuesX"][0]
        for j in range(len(d["categories"])):
            axs[i].plot(
                x,
                d["valuesY"][j],
                label=d["categories"][j],
                color=d["colors"][j],
                marker='o',        # shape: 'o', 's', '^', 'D', etc.
                markersize=5,      # size of the markers
                markeredgecolor='black',  # outline (optional)
                linewidth=1.5      # thickness of the connecting line
            )

        axs[i].set_title(d["title"], fontsize=12, fontweight='bold')
        axs[i].set_ylabel(yLable)
        axs[i].grid(axis='y', linestyle='--', alpha=0.5)
        axs[i].legend()

def barGraphs(fullDict,baseDir,imgTitle,yLable,config,titles):
    graphsNo = len(config)
    data = baseDicts(graphsNo,["categories","values","title","colors"])
    addTitles(data,titles)

    for i in range(len(config)):
        for bar in config[i]:
            addBar(fullDict,bar[0],bar[1],data,i)

    addColors(data,"values")

    fig, axs = plt.subplots(graphsNo, 1, figsize=(10, 5*graphsNo))

    axs = [axs] if graphsNo == 1 else axs.flatten()

    plotBarData(data,axs,yLable)

    fig.suptitle(imgTitle, fontsize=16, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(f"{baseDir}/{imgTitle}.png", dpi=300, bbox_inches="tight")

def lineGraphs(fullDict,baseDir,imgTitle,yLable,config,titles):
    graphsNo = len(config)
    data = baseDicts(graphsNo,["categories","valuesX","valuesY","title","colors"])
    addTitles(data,titles)

    for i in range(len(config)):
        for line in config[i]:
            addLine(fullDict,line[0],line[1],data,i)
    
    addColors(data,"valuesY")

    fig, axs = plt.subplots(graphsNo, 1, figsize=(10, 5*graphsNo))

    axs = [axs] if graphsNo == 1 else axs.flatten()

    plotLineData(data,axs,yLable)

    fig.suptitle(imgTitle, fontsize=16, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(f"{baseDir}/{imgTitle}.png", dpi=300, bbox_inches="tight")


def getDict(filePath):
    with open(filePath) as json_file:
        data = json.load(json_file)

    return data

def plotGroup(saveFilesBasePath,groupName,groupPath):
    fullDict = getDict(groupPath)
    dirPath = f"{saveFilesBasePath}/{groupName}_Plots"

    if os.path.exists(dirPath):
        shutil.rmtree(dirPath)

    os.makedirs(dirPath,exist_ok=True)

    barGraphs(fullDict,dirPath,f"{groupName} Percentage Statistics","Percentage",percentageConfig,percentageTitles)
    barGraphs(fullDict,dirPath,f"{groupName} Exon Size Statistics","Nucleotide Bases",exonSizeConfig,exonSizeTitles)
    barGraphs(fullDict,dirPath,f"{groupName} Intron Size Statistics","Nucleotide Bases",intronSizeConfig,intronSizeTitles)
    barGraphs(fullDict,dirPath,f"{groupName} Exon Quantity Statistics","Number of Exons",avgExonsInGenesConfig,avgExonsInGenesTitle)
    lineGraphs(fullDict,dirPath,f"{groupName} Group Percentage Statistics","Percentage",groupPercentageConfig,groupPercentageTitles)
    lineGraphs(fullDict,dirPath,f"{groupName} Group Size Statistics","Average Size Percentil",groupSizeConfig,groupSizeTitles)


def plot(saveFilesBasePath):
    plotGroup(saveFilesBasePath,"Precision",f"{saveFilesBasePath}/finalJsons/precisionStatistics.json")
    plotGroup(saveFilesBasePath,"Recall",f"{saveFilesBasePath}/finalJsons/recallStatistics.json")
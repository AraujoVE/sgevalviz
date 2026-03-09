import matplotlib.pyplot as plt
import json
import os
import shutil

percentageConfig = [
    [
        ("nucleotides_predicted_percentage_forward","Forward Strand"),
        ("nucleotides_predicted_percentage_reverse","Reverse Strand"),
        ("nucleotides_predicted_percentage","Total")
    ],
    [
        ("genes_predicted_percentage_forward","Forward Strand"),
        ("genes_predicted_percentage_reverse","Reverse Strand"),
        ("genes_predicted_percentage","Total")
    ],
    [
        ("start_codons_predicted_percentage_forward","Forward Strand"),
        ("start_codons_predicted_percentage_reverse","Reverse Strand"),
        ("start_codons_predicted_percentage","Total")
    ],
    [
        ("stop_codons_predicted_percentage_forward","Forward Strand"),
        ("stop_codons_predicted_percentage_reverse","Reverse Strand"),
        ("stop_codons_predicted_percentage","Total")
    ],
    [
        ("genes_not_ignored_percentage_forward", "Forward Strand"),
        ("genes_not_ignored_percentage_reverse", "Reverse Strand"),
        ("genes_not_ignored_percentage", "Total")
    ],
    [
        ("exons_predicted_percentage_forward","Forward Strand"),
        ("exons_predicted_percentage_reverse","Reverse Strand"),
        ("exons_predicted_percentage","Total")
    ],
    [
        ("first_exons_only_predicted_percentage_forward","Forward Strand"),
        ("first_exons_only_predicted_percentage_reverse","Reverse Strand"),
        ("first_exons_only_predicted_percentage","Total")
    ],
    [
        ("last_exons_only_predicted_percentage_forward","Forward Strand"),
        ("last_exons_only_predicted_percentage_reverse","Reverse Strand"),
        ("last_exons_only_predicted_percentage","Total")
    ],
    [
        ("single_exons_predicted_percentage_forward","Forward Strand"),
        ("single_exons_predicted_percentage_reverse","Reverse Strand"),
        ("single_exons_predicted_percentage","Total"),
    ],
    [
        ("introns_predicted_percentage_forward","Forward Strand"),
        ("introns_predicted_percentage_reverse","Reverse Strand"),
        ("introns_predicted_percentage","Total")
    ],
    [
        ("donors_predicted_percentage_forward","Forward Strand"),
        ("donors_predicted_percentage_reverse","Reverse Strand"),
        ("donors_predicted_percentage","Total")
    ],
    [
        ("acceptors_predicted_percentage_forward","Forward Strand"),
        ("acceptors_predicted_percentage_reverse","Reverse Strand"),
        ("acceptors_predicted_percentage","Total")
    ]
]
percentageTitles = ("Nucleotides", "Genes", "Start Codon", "Stop Codon", "Genes Not Ignored", "General Exons", "First Exons", "Last Exons", "Single Exons", "Introns", "Donors", "Acceptors")

exonSizeConfig = [
    [
        ("exons_total_size_avg_forward","All, Forward Strand"),
        ("exons_total_size_avg_reverse","All, Reverse Strand"),
        ("exons_total_size_avg","All, Total"),
        ("exons_predicted_size_avg_forward","Predicted, Forward Strand"),
        ("exons_predicted_size_avg_reverse","Predicted, Reverse Strand"),
        ("exons_predicted_size_avg","Predicted, Total"),
        ("exons_unpredicted_size_avg_forward","Unpredicted, Forward Strand"),
        ("exons_unpredicted_size_avg_reverse","Unpredicted, Reverse Strand"),
        ("exons_unpredicted_size_avg","Unpredicted, Total")
    ],
    [
        ("first_exons_only_total_size_avg_forward","All, Forward Strand"),
        ("first_exons_only_total_size_avg_reverse","All, Reverse Strand"),
        ("first_exons_only_total_size_avg","All, Total"),
        ("first_exons_only_predicted_size_avg_forward","Predicted, Forward Strand"),
        ("first_exons_only_predicted_size_avg_reverse","Predicted, Reverse Strand"),
        ("first_exons_only_predicted_size_avg","Predicted, Total"),
        ("first_exons_only_unpredicted_size_avg_forward","Unpredicted, Forward Strand"),
        ("first_exons_only_unpredicted_size_avg_reverse","Unpredicted, Reverse Strand"),
        ("first_exons_only_unpredicted_size_avg","Unpredicted, Total")
    ],
    [
        ("last_exons_only_total_size_avg_forward","All, Forward Strand"),
        ("last_exons_only_total_size_avg_reverse","All, Reverse Strand"),
        ("last_exons_only_total_size_avg","All, Total"),
        ("last_exons_only_predicted_size_avg_forward","Predicted, Forward Strand"),
        ("last_exons_only_predicted_size_avg_reverse","Predicted, Reverse Strand"),
        ("last_exons_only_predicted_size_avg","Predicted, Total"),
        ("last_exons_only_unpredicted_size_avg_forward","Unpredicted, Forward Strand"),
        ("last_exons_only_unpredicted_size_avg_reverse","Unpredicted, Reverse Strand"),
        ("last_exons_only_unpredicted_size_avg","Unpredicted, Total")
    ],
    [
        ("single_exons_total_size_avg_forward","All, Forward Strand"),
        ("single_exons_total_size_avg_reverse","All, Reverse Strand"),
        ("single_exons_total_size_avg","All, Total"),
        ("single_exons_predicted_size_avg_forward","Predicted, Forward Strand"),
        ("single_exons_predicted_size_avg_reverse","Predicted, Reverse Strand"),
        ("single_exons_predicted_size_avg","Predicted, Total"),
        ("single_exons_unpredicted_size_avg_forward","Unpredicted, Forward Strand"),
        ("single_exons_unpredicted_size_avg_reverse","Unpredicted, Reverse Strand"),
        ("single_exons_unpredicted_size_avg","Unpredicted, Total")
    ]
]
exonSizeTitles = ["General Exon Size","General First Exon Size","General Last Exon Size","General Single Exon Size"]

intronSizeConfig = [
    [
        ("introns_total_size_avg_forward","All, Forward Strand"),
        ("introns_total_size_avg_reverse","All, Reverse Strand"),
        ("introns_total_size_avg","All, Total"),
        ("introns_predicted_size_avg_forward","Predicted, Forward Strand"),
        ("introns_predicted_size_avg_reverse","Predicted, Reverse Strand"),
        ("introns_predicted_size_avg","Predicted, Total"),
        ("introns_unpredicted_size_avg_forward","Unpredicted, Forward Strand"),
        ("introns_unpredicted_size_avg_reverse","Unpredicted, Reverse Strand"),
        ("introns_unpredicted_size_avg","Unpredicted, Total")
    ]
]
intronSizeTitles = ["General Intron"]


avgExonsInGenesConfig = [
    [
        ("exons_in_predicted_genes_predicted_transcripts_forward","Forward Strand"),
        ("exons_in_predicted_genes_predicted_transcripts_reverse","Reverse Strand"),
        ("exons_in_predicted_genes_predicted_transcripts","Total")
    ],
    [
        ("exons_in_predicted_genes_unpredicted_transcripts_forward","Forward Strand"),
        ("exons_in_predicted_genes_unpredicted_transcripts_reverse","Reverse Strand"),
        ("exons_in_predicted_genes_unpredicted_transcripts","Total")
    ],
    [
        ("exons_in_unpredicted_genes_unpredicted_transcripts_forward","Forward Strand"),
        ("exons_in_unpredicted_genes_unpredicted_transcripts_reverse","Reverse Strand"),
        ("exons_in_unpredicted_genes_unpredicted_transcripts","Total")
    ]
]
avgExonsInGenesTitle = ["Average Number of Exons in Predicted Genes, Predicted Transcripts","Average Number of Exons in Predicted Genes, Unpredicted Transcripts","Average Number of Exons in Unredicted Genes"]

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

def getDict(filePath):
    with open(filePath) as json_file:
        data = json.load(json_file)

    return data

def plotGroup(saveFilesBasePath,groupName,groupPath):
    fullDict = getDict(groupPath)
    dirPath = f"{saveFilesBasePath}{groupName}_Plots"

    if os.path.exists(dirPath):
        shutil.rmtree(dirPath)

    os.makedirs(dirPath,exist_ok=True)

    barGraphs(fullDict,dirPath,f"{groupName} Percentage Statistics","Percentage",percentageConfig,percentageTitles)
    #barGraphs(fullDict,dirPath,f"{groupName} Exon Size Statistics","Nucleotide Bases",exonSizeConfig,exonSizeTitles)
    #barGraphs(fullDict,dirPath,f"{groupName} Intron Size Statistics","Nucleotide Bases",intronSizeConfig,intronSizeTitles)
    barGraphs(fullDict,dirPath,f"{groupName} Exon Quantity Statistics","Number of Exons",avgExonsInGenesConfig,avgExonsInGenesTitle)


def plot(saveFilesBasePath):
    plotGroup(saveFilesBasePath,"Precision",f"{saveFilesBasePath}finalJsons/precisionStatistics.json")
    plotGroup(saveFilesBasePath,"Recall",f"{saveFilesBasePath}finalJsons/recallStatistics.json")
#!/usr/bin/python3

import glob
import subprocess
import argparse
from pathlib import Path

import json
import os

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, )
parser.add_argument("-o", "--output", type=str, default="")
p = parser.parse_args()

infolder = Path(p.input).resolve()
database = Path(p.database).resolve()

# Create an output directory
outputDir = p.output
inputFile = p.input
os.mkdir(outputDir)


# Check input
def checkin(inputfolder):

    for fastfolder in inputfolder:
        print("FASTQ-folder:", fastfolder)

    return inputfolder


def run(inputfolders, outputfolder):
    for fastqfolder in inputfolders:
        foldername = fastqfolder["fastqFolder"]
        sid = fastqfolder["ID"]
        outputfolder = outputfolder + str(sid)
        if fastqfolder["gzip"]:
            filetype = ".fastq.gz"
        else:
            filetype = ".fastq"
        phage = fastqfolder["phageGenome"]
        genomesize = fastqfolder["genomesize"]
        pipeline = glob.glob("pipeline.sh")[0]
        if phage:
            kingdom = "Virus"
        else:
            kingdom = "Bacteria"
        subprocess.run(["bash", pipeline, "-i", foldername, "-o", outputfolder, "-g", genomesize, "-c", "50",
                        "-t", "1000", "-q", "12", "-k", kingdom, "-f", filetype], text=True)


# main
folderlist = checkin(assemblyfolders)
run(folderlist, outputDir)

outputJson = {
    'resultsfolder': outputDir
}
with open('/batchx/output/output.json', 'w+') as json_file:
    json.dump(outputJson, json_file)


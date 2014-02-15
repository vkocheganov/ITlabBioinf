class probeClass:
	pass
def ExtractProbeInfo(fileName):
	probes = []
	with open(fileName)	as csvfile:
		fileReader = csv.reader(csvfile, delimiter='\t')
		lineInd = 0
		for	line in fileReader:
			if lineInd < 2:
				lineInd += 1
				continue
			wordInd = 0
			myProbe = probeClass()
			myProbe.ind = 0
			for	word in line:
				if ((word == 'NA') or (word == '') or (word.find(';') != -1 )) and ( wordInd == 0 or wordInd == 2 or wordInd == 4):
					myProbe.ind = -1
				if wordInd == 0:
					myProbe.cgName = word
				elif wordInd == 2:
					myProbe.geneName = word
				elif wordInd == 4:
					myProbe.geneCoord = word
				wordInd += 1
			probes.append(myProbe)
			lineInd += 1

	return probes

def transpose(fileName, outFileName):
	samples = []
	with open(fileName) as csvfile:
		fileReader = csv.reader(csvfile, delimiter=' ')
		for line in fileReader:
			samples.append(line)
#			print(len(samples))
	fileWriter = open(outFileName,'w')
	for wordCount in range(len(samples[0]) - 1):
		for lineCount in range(len(samples) - 1):
			fileWriter.write(samples[lineCount][wordCount] + " ")
		fileWriter.write(samples[len(samples) - 1][wordCount] + "\n")
#		print(wordCount)
	fileWriter.close()
#import csv
#transpose("c:/Users/Victor/Downloads/samples.csv","mysamples.csv")
#exit(0)
def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False

def ExtractProbeAvailInfo(fileName):
	avail = []
	with open(fileName)	as csvfile:
		fileReader = csv.reader(csvfile, delimiter='\t')
		lineInd = 0
		for	line in fileReader:
			if lineInd < 2:
				lineInd += 1
				continue
			wordInd = 0
			avail.append(int(line[1].find('.') != -1))

	return avail

def mergeAvails(avail0, avail1):
	count = 0
	for ind in range(len(avail1)):
		avail0[ind] = int((avail0[ind] == 1) and (avail1[ind] == 1))


def mergeProbeWithAvail(probes, avail):
	for ind in range(len(avail)):
		if (avail[ind] == 0 or probes[ind].ind == -1):
			probes[ind].ind = -1

def WriteProbeInfo(probes, outFileName):
	with open(outFileName,'w+')	as outCsvFile:
		fileWriter = csv.writer(outCsvFile, delimiter=' ', lineterminator="\n")
		for probe in probes:
			fileWriter.writerow([probe.cgName, probe.geneName, probe.geneCoord,	probe.ind])

def ReadProbeInfo(fileName):
	probes = []
	availProbes = 0
	with open(fileName)	as csvfile:
		fileReader = csv.reader(csvfile, delimiter=' ')
		lineInd = 0
		for	line in fileReader:
			myProbe = probeClass()
			myProbe.cgName = line[0]
			myProbe.geneName = line[1]
			myProbe.geneCoord = line[2]
			myProbe.ind = line[3]
			if line[3] != "-1":
				availProbes += 1
			probes.append(myProbe)
	return (availProbes, probes)

def ProcessSample(fileWriter, probes, availProbes, inputFileName):
	sample = [0 for x in range(availProbes)]
	ind = 0
	with open(inputFileName) as csvfile:
		fileReader = csv.reader(csvfile, delimiter='\t')
		for line in fileReader:
			if ind < 2:
				ind += 1
				continue
			sampleInd = int(probes[ind - 2].ind)
			if (sampleInd != -1):
				sample[sampleInd] = float(line[1])
			ind += 1
	for probe in sample:
		fileWriter.write(str(probe) + " ")
	fileWriter.write("\n")

def isHealphy(fileName):
#	print(fileName.split("-")[5][0:2])
	if int(fileName.split("-")[5][0:2]) < 10:
		return 0
	else:
		return 1

def CreateSamples(folderToRead, firstFile, outFileNameSamples, outFileNameResponses):
	outFileName	= "c:/Users/Victor/Downloads/newProbeInfo.csv"
	probes = ExtractProbeInfo(firstFile)
	avail0 = ExtractProbeAvailInfo(firstFile)
	avail1 = []

	print("processing folder for avail")
	for dirpath, dnames, fnames in os.walk(folderToRead):
		for	f in fnames:
			if f[0:3] == "jhu":
				print("process ", os.path.join(dirpath, f))
				avail1 = ExtractProbeAvailInfo(os.path.join(dirpath, f))
				mergeAvails(avail0,avail1)
	print("merging probes with avail")
	mergeProbeWithAvail(probes,avail0)
	print("writing probes to file")
	WriteProbeInfo(probes, outFileName)

	print("reading probes from file")
	(availProbes, probes) =	ReadProbeInfo(outFileName)
	mysort = sorted(range(len(probes)), key = lambda y: probes[y].geneName)
	count = 0
	for i in mysort:
		if (probes[i].ind != "-1"):
			probes[i].ind =	count
			count += 1
	fileWriterSamples = open(outFileNameSamples,"w")
	fileWriterResponses = open(outFileNameResponses,"w")
	print("forming samples file")
	for dirpath, dnames, fnames in os.walk(folderToRead):
		for	f in fnames:
			if f[0:3] == "jhu":
				print("process ", os.path.join(dirpath, f))
				ProcessSample(fileWriterSamples, probes, availProbes, os.path.join(dirpath, f))
				fileWriterResponses.write(str(1-isHealphy(f)) + " ")
	fileWriterSamples.close()

import csv
import array
import time
import os
#from itertools import groupby
#probefileName = "c:/Users/Victor/Documents/betaData/extracted/jhu-usc.edu_BLCA.HumanMethylation450.Level_3.1.12.0/jhu-usc.edu_BLCA.HumanMethylation450.1.lvl-3.TCGA-AV-A03D-20A-01D-A10W-05.txt "
#probefileName = "c:/Users/Victor/Documents/betaData/KIRC/extracted/jhu-usc.edu_KIRC.HumanMethylation450.Level_3.1.9.0/jhu-usc.edu_KIRC.HumanMethylation450.1.lvl-3.TCGA-B0-4813-01A-01D-1275-05.txt"
probefileName = "c:/Users/Victor/Documents/betaData/BLCA/Extracted/jhu-usc.edu_BLCA.HumanMethylation450.Level_3.1.12.0/jhu-usc.edu_BLCA.HumanMethylation450.1.lvl-3.TCGA-AV-A03D-20A-01D-A10W-05.txt"
outFileNameSamples = "c:/Users/Victor/Downloads/samples.csv"
outFileNameResponses = "c:/Users/Victor/Downloads/responses.csv"
start = time.time()
#CreateSamples("c:/Users/Victor/Documents/betaData/extracted/",probefileName, outFileNameSamples, outFileNameResponses)
CreateSamples("c:/Users/Victor/Documents/betaData/BLCA/Extracted/",probefileName, outFileNameSamples, outFileNameResponses)
transpose(outFileNameSamples,outFileNameSamples)
print("time spent:", time.time()-start," secs")




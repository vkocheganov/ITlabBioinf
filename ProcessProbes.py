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
			notAvail = 0
			for	word in line:
				if ((word == 'NA') or (word == '') or (word.find(';') != -1 )) and ( wordInd == 0 or wordInd == 2 or wordInd == 4):
					notAvail = 1
				wordInd += 1
			wordInd = 0
			if notAvail == 0:
				myProbe = probeClass()
				for	word in line:
					if wordInd == 0:
						myProbe.cgName = word
					elif wordInd == 2:
						myProbe.geneName = word
					elif wordInd == 4:
						myProbe.geneCoord = word
					wordInd += 1
				myProbe.ind = 0
				probes.append(myProbe)
			lineInd += 1
	return probes
def WriteProbeInfo(probes, indeces, outFileName):
	with open(outFileName,'w+')	as outCsvFile:
		fileWriter = csv.writer(outCsvFile, delimiter=' ')
#		probeInd = 0
		for probe in probes:
			#fileWriter.writerow([probe.cgName, probe.geneName, probe.geneCoord, indeces[probeInd]])
			fileWriter.writerow([probe.cgName, probe.geneName, probe.geneCoord, probe.ind])
#			probeInd += 1
			
import csv;
import array;
import time;

#from itertools import groupby
probefileName = "c:/Users/Victor/Downloads/jhu-usc.edu_BLCA.HumanMethylation450.Level_3.1.12.0/jhu-usc.edu_BLCA.HumanMethylation450.1.lvl-3.TCGA-AV-A03D-20A-01D-A10W-05.txt"
outFileName = "c:/Users/Victor/Downloads/newProbeInfo.csv"

start = time.time()

probes = ExtractProbeInfo(probefileName)
print("probes count = ", len(probes))

mysort = sorted(range(len(probes)), key = lambda y: probes[y].geneName)
count = 0

for i in mysort:
	probes[i].ind = count
	count += 1

WriteProbeInfo(probes, [],outFileName)
print("time spent:", time.time()-start," secs")
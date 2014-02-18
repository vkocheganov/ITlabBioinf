import csv
import array
import time
import os
import sys

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
                                # info of interest is stored in 0, 2 and 4 words:
                                # 0 - CpG site name
                                # 2 - corresponding gene name
                                # 4 - corresponding geneCoord
                                # Skip info (set ind to -1) if:
                                # if info is unavailable for CpG site
                                # or info is equal to ''
                                # or there are several Genese correspond to this site (they are ';' separated)
                                # 
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

def ExtractProbeAvailInfo(fileName, size):
        # Check for this particular file what CpG sites have beta value available
        # This info is file-specific
        avail = []
        with open(fileName)	as csvfile:
                fileReader = csv.reader(csvfile, delimiter='\t')
                next(fileReader,None)
                next(fileReader,None)
                for line in fileReader:
                        avail.append(int(line[1].find('.') != -1))
	return avail

def mergeAvails(avail0, avail1):
        # element-wise operation:
        # avail0 = avail0 && avail1 
        print "merging\n"
	count = 0
	for ind in range(len(avail1)):
		avail0[ind] = int((avail0[ind] == 1) and (avail1[ind] == 1))

def mergeProbeWithAvail(probes, avail):
	for ind in range(len(avail)):
		if (avail[ind] == 0 or probes[ind].ind == -1):
			probes[ind].ind = -1

def WriteProbeInfo(probes, outFileName):
        # Be aware about ReadProbeInfo
	with open(outFileName,'w+')	as outCsvFile:
		fileWriter = csv.writer(outCsvFile, delimiter=' ', lineterminator="\n")
		for probe in probes:
			fileWriter.writerow([probe.cgName, probe.geneName, probe.geneCoord,	probe.ind])

def ReadProbeInfo(fileName):
        # Be aware about WriteProbeInfo
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

def ProcessSample(probes, availProbes, inputFileName):
        print "Process samples\n"
        start = time.time()
	sample = [0 for x in range(availProbes)]
        ind = 0
	with open(inputFileName) as csvfile:
		fileReader = csv.reader(csvfile, delimiter='\t')
                next(fileReader)
                next(fileReader)
		for line in fileReader:
			sampleInd = int(probes[ind].ind)
			if (sampleInd != -1):
				sample[sampleInd] = float(line[1])
			ind += 1
        
        print("time spent:", time.time()-start," secs")
        return sample

def isHealphy(fileName):
	if int(fileName.split("-")[5][0:2]) < 10:
		return 0
	else:
		return 1

def CreateSamples(folderToRead, firstFile, outFileNameSamples, outFileNameResponses, outFileName):
        print "Creating samples"
        print "Get only CpG sites with sufficient accompanying info\n:"
	probes = ExtractProbeInfo(firstFile)
        
        # Very similar procedure, but analize CpG sites in sense of available beta value;
        # since this "availability" differs from file to file new function was added:
	avail0 = ExtractProbeAvailInfo(firstFile,0)
	avail1 = []
	print "processing folder to get available CpG sites"
        filecount = 0
	for dirpath, dnames, fnames in os.walk(folderToRead):
		for	f in fnames:
			if (f[0:3] == "jhu") and (os.path.join(dirpath, f) != firstFile):
                                filecount +=1
				print "process " ,filecount, " ", os.path.join(dirpath, f)
				avail1 = ExtractProbeAvailInfo(os.path.join(dirpath, f), len(avail0))
				mergeAvails(avail0,avail1)
                
	print("merging probes with avail")
        start = time.time()
	mergeProbeWithAvail(probes,avail0)
        print("time spent:", time.time()-start," secs")
	print("writing probes to file")
        start = time.time()
	WriteProbeInfo(probes, outFileName)
        print("time spent:", time.time()-start," secs")

	print("reading probes from file")
        start = time.time()
	(availProbes, probes) =	ReadProbeInfo(outFileName)
        print("time spent:", time.time()-start," secs")
        print "sorting\n"
        start = time.time()
	mysort = sorted(range(len(probes)), key = lambda y: int(probes[y].geneCoord))
	count = 0
	for i in mysort:
		if (probes[i].ind != "-1"):
			probes[i].ind =	count
			count += 1
	WriteProbeInfo(probes, outFileName)
        print("time spent:", time.time()-start," secs")
	fileWriterSamples = open(outFileNameSamples,"w")
	fileWriterResponses = open(outFileNameResponses,"w")
	print("forming samples file")
        start = time.time()

	samples = []
        filecount = 0
        for dirpath, dnames, fnames in os.walk(folderToRead):
		for	f in fnames:
			if (f[0:3] == "jhu"):
				print("process ", os.path.join(dirpath, f))
                                filecount +=1
				samples.append(ProcessSample(probes, availProbes, os.path.join(dirpath, f)))
				fileWriterResponses.write(str(1-isHealphy(f)) + " ")
	for wordCount in range(len(samples[0])):
		for lineCount in range(len(samples) - 1):
			fileWriterSamples.write(str(samples[lineCount][wordCount]) + " ")
		fileWriterSamples.write(str(samples[len(samples) - 1][wordCount]) + "\n")
	fileWriterSamples.close()
        
        print("time spent:", time.time()-start," secs")

def ProcessCancer(folderToProcess, workFolder,cancer):
        print "############# Process cancer " + cancer
        start = time.time()
        cancerFolder=folderToProcess + "/" + cancer + "/Extracted/"
        firstProbefileName=""
        for dirpath, dnames, fnames in os.walk(cancerFolder):
                for f in fnames:
                        if f[0:3] == "jhu":
                                firstProbefileName=os.path.join(dirpath, f)
                                print "first file: "+firstProbefileName
                                break
                        if (firstProbefileName != ""):
                                break
                if (firstProbefileName != ""):
                        break
        outFileNameSamples = workFolder+"/" + cancer + ".sample.csv"
        outFileNameResponses = workFolder+"/" + cancer + ".resp.csv"
        print "outFileNameSamples: " +outFileNameSamples
        print "outFileNameResponses: " + outFileNameResponses
	outFileName = workFolder+"/"+ cancer + "newProbeInfo.csv"
        CreateSamples(cancerFolder,firstProbefileName, outFileNameSamples, outFileNameResponses,outFileName)
        print("time spent:", time.time()-start," secs")

# Main part
if (len(sys.argv) < 3):
        print("Please provide with directory to process and to work with")
        exit(1)

folderToProcess=sys.argv[1]
print "folderToProcess: " + folderToProcess
workFolder=sys.argv[2]

cancers = []
cancers.append("BLCA");
cancers.append("READ");
cancers.append("KIRP");
cancers.append("LIHC");
cancers.append("PRAD");
cancers.append("LUSC");
cancers.append("COAD");
cancers.append("LUAD");
cancers.append("HNSC");
# cancers.append("THCA");
# cancers.append("UCEC");
# cancers.append("KIRC");
# cancers.append("BRCA");
for cancer in cancers:
        ProcessCancer(folderToProcess,workFolder,cancer)
        
exit(0)

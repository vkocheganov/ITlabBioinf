from numpy import *
from numpy.linalg import *
from scipy.stats import norm

trainRatio=0.5
smoothingEps=0.0000000000000000000000000001
ExperimentNumber=1
folderName="/home/victor/Documents/beta_data/ITlabBioinf/temp/"
#folderName="G:/betta_data/temp"
cancerName="READ"

def MakeExperiment(hSamples, hTrainSize, hTestSize,cSamples, cTrainSize, cTestSize):
    hrandPerm=random.permutation(len(hSamples[0]))
    #hTrainInd=hrandPerm[range(hTrainSize)]
    hTrainInd=hrandPerm[array(range(hTrainSize))]
    hTestInd=hrandPerm[array([hTrainSize]) + range(hTestSize)]
    hTrainSamples=hSamples[:,hTrainInd]
    hTestSamples=hSamples[:,hTestInd]
    hMeans=mean(hTrainSamples,1)
    hsds=std(hTrainSamples,1)
    hCorrs=zeros(len(hMeans)-1)
    for ind in range(len(hMeans)-1):
        hCorrs[ind] = corrcoef(hTrainSamples[ind,:],hTrainSamples[ind+1,:])[0,1]
    crandPerm=random.permutation(len(cSamples[0]))
    # cTrainInd=crandPerm[range(cTrainSize)]
    cTrainInd=crandPerm[array(range(cTrainSize))]
    cTestInd=crandPerm[array([cTrainSize]) + range(cTestSize)]
    cTrainSamples=cSamples[:,cTrainInd]
    cTestSamples=cSamples[:,cTestInd]
    cMeans=mean(cTrainSamples,1)
    csds=std(cTrainSamples,1)
    cCorrs=zeros(len(cMeans)-1)
    for ind in range(len(cMeans)-1):
        cCorrs[ind] = corrcoef(cTrainSamples[ind,:],cTrainSamples[ind+1,:])[0,1]
    herrs = zeros(hTestSize,dtype=bool)
    for ind in range(hTestSize):
        hCondMeans = hMeans[1:len(hMeans)] + divide( multiply( multiply( hsds[1:len(hsds)],hCorrs ),hTestSamples[0:(len(hMeans)-1),ind]-hMeans[0:(len(hMeans)-1)] ),hsds[0:(len(hsds)-1)] )
        hCondSds = multiply(sqrt(ones(len(hMeans)-1) - square(hCorrs)),hsds[1:len(hMeans)])
        prob0=sum(log(norm.pdf(divide(hTestSamples[1:len(hMeans),ind] - hCondMeans ,hCondSds))+smoothingEps)) + log(float(hTrainSize)/(hTrainSize+cTrainSize)) + log(norm.pdf((hTestSamples[0,ind] - hMeans[0])/hsds[0]) + smoothingEps)
        cCondMeans = cMeans[1:len(cMeans)] + divide( multiply( multiply( csds[1:len(hsds)],cCorrs ),hTestSamples[0:(len(cMeans)-1),ind]-cMeans[0:(len(cMeans)-1)] ),csds[0:(len(csds)-1)] )
        cCondSds = multiply(sqrt(ones(len(cMeans)-1) - square(cCorrs)),csds[1:len(cMeans)])
        prob1=sum(log(norm.pdf(divide(hTestSamples[1:len(cMeans),ind] - cCondMeans,cCondSds))+smoothingEps)) + log(float(cTrainSize)/(hTrainSize+cTrainSize))+ log(norm.pdf((hTestSamples[0,ind] - cMeans[0])/csds[0]) + smoothingEps)
        # prob0=sum(log(norm.pdf(divide(hTestSamples[:,ind] - hMeans,hsds))+smoothingEps)) + hTrainSize/(hTrainSize+cTrainSize)
        # prob1=sum(log(norm.pdf(divide(hTestSamples[:,ind] - cMeans,csds))+smoothingEps)) + cTrainSize/(hTrainSize+cTrainSize)
        # print( "prob0 = ",prob0," prob1 = ", prob1)
        herrs[ind] = prob0<prob1
    cerrs = zeros(cTestSize,dtype=bool)
    for ind in range(cTestSize):
        hCondMeans = hMeans[1:len(hMeans)] + divide( multiply( multiply( hsds[1:len(hsds)],hCorrs ),cTestSamples[0:(len(hMeans)-1),ind]-hMeans[0:(len(hMeans)-1)] ),hsds[0:(len(hsds)-1)] )
        hCondSds = multiply(sqrt(ones(len(hMeans)-1) - square(hCorrs)),hsds[1:len(hMeans)])
        prob0=sum(log(norm.pdf(divide(cTestSamples[1:len(hMeans),ind] - hCondMeans ,hCondSds))+smoothingEps)) + log(float(hTrainSize)/(hTrainSize+cTrainSize))+log(norm.pdf((cTestSamples[0,ind] - hMeans[0])/hsds[0]) + smoothingEps)
        cCondMeans = cMeans[1:len(cMeans)] + divide( multiply( multiply( csds[1:len(hsds)],cCorrs ),cTestSamples[0:(len(cMeans)-1),ind]-cMeans[0:(len(cMeans)-1)] ),csds[0:(len(csds)-1)] )
        cCondSds = multiply(sqrt(ones(len(cMeans)-1) - square(cCorrs)),csds[1:len(cMeans)])
        prob1=sum(log(norm.pdf(divide(cTestSamples[1:len(cMeans),ind] - cCondMeans,cCondSds))+smoothingEps)) + log(float(cTrainSize)/(hTrainSize+cTrainSize))+ log(norm.pdf((cTestSamples[0,ind] - cMeans[0])/csds[0]) + smoothingEps)
        # prob0=sum(log(norm.pdf(divide(cTestSamples[:,ind] - hMeans,hsds))+smoothingEps)) + hTrainSize/(hTrainSize+cTrainSize)
        # prob1=sum(log(norm.pdf(divide(cTestSamples[:,ind] - cMeans,csds))+smoothingEps)) + cTrainSize/(hTrainSize+cTrainSize)
        # print( "prob0 = ",prob0," prob1 = ", prob1)
        cerrs[ind] = prob0>prob1
    return sum(herrs),sum(cerrs)

def ProcessCancer(folderName,cancerName,filewrite):
    print("Processing ", cancerName)
    sampleFileName=folderName+"/"+cancerName+".sample.csv"
    respFileName=folderName+"/"+cancerName+".resp.csv"
    samples= loadtxt(sampleFileName, dtype=float)
    responses= loadtxt(respFileName, dtype=bool)
    hSamples=samples[:,~responses]
    hTrainSize=int(len(hSamples[0])*trainRatio)
    hTestSize=len(hSamples[0])-hTrainSize
    cSamples=samples[:,responses]
    cTrainSize=int(len(cSamples[0])*trainRatio)
    cTestSize=len(cSamples[0])-cTrainSize
    sumHErr=0
    sumCErr=0
    for ind in range(ExperimentNumber):
        herr,cerr=MakeExperiment(hSamples, hTrainSize, hTestSize,cSamples, cTrainSize, cTestSize)
        print("ind = ",ind," herr = ",herr,"/",hTestSize, " cerr = ",cerr,"/",cTestSize)
        sumHErr+=herr
        sumCErr+=cerr
    print("Total errs: \n")
    herr = float(sumHErr)/(hTestSize*ExperimentNumber)
    cerr = float(sumCErr)/(cTestSize*ExperimentNumber)
    print("Healphy: ",sumHErr,"/",hTestSize*ExperimentNumber, " (",herr,")\n")
    print("Cancer: ",sumCErr,"/",cTestSize*ExperimentNumber, " (",cerr,")\n")
    filewrite.write(cancerName + " " + str(ExperimentNumber) + " " + str(len(hSamples[0])) + " "+ str(len(cSamples[0])) + " " + str(herr) + " " + str(cerr))
    
cancers = []
# cancers.append("BLCA");
cancers.append("READ");
# cancers.append("KIRP");
# cancers.append("LIHC");
# cancers.append("PRAD");
# cancers.append("LUSC");
# cancers.append("COAD");
# cancers.append("LUAD");
# cancers.append("HNSC");
# cancers.append("THCA");
# cancers.append("UCEC");
# cancers.append("KIRC");
# cancers.append("BRCA");
filewrite = open("MarkovReport.out","w")
filewrite.write("cancerName " + "NExp " + "HealphyNum " +  "CancerNum " + "HealphyErr "+ "CancerErr\n" )
for cancer in cancers:
    ProcessCancer(folderName,cancer,filewrite)
filewrite.close()

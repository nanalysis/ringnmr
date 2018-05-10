import org.comdnmr.cpmgfit2.calc.CalcRDisp as CalcRDisp
import org.comdnmr.cpmgfit2.calc.DataIO as DataIO
import org.comdnmr.cpmgfit2.calc.ResidueFitter as ResidueFitter
import math
import sys
import os.path
import runpy
#from optparse import OptionParser
import argparse

calcR = CalcRDisp()

def findDuplicates(vcpmgs):
    dupDict = {}
    for i,vcpmg in enumerate(vcpmgs):
        if not vcpmg in dupDict:
            dupDict[vcpmg] = []
        dupDict[vcpmg].append(i)
   
    dups = []
    for i,vcpmg in enumerate(dupDict):
        if len(dupDict[vcpmg]) > 1:
            dups.append(dupDict[vcpmg])
    return dups

def getSDevFromPeakDups(duplicates,usedValues,useY):
    sdev = 0.0
    sumDelta2 = 0.0
    nDup = 0
    for dup in duplicates:
        (dup1,dup2) = dup
        udup1 = usedValues[dup1]
        udup2 = usedValues[dup2]
        if udup1 != -1 and udup2 != -1:
            v1 = useY[udup1]
            v2 = useY[udup2]
            delta = v1-v2
            sumDelta2 += delta*delta
            nDup += 1

    if nDup > 0:
        sdev = math.sqrt(sumDelta2/(2*nDup))
    return sdev

def dumpData(func,x,y,err,fields,fieldValues,pars):
    func.setXY(x,y)
    func.setErr(err)
    func.setFieldValues(fieldValues)
    func.setFields(fields)
    func.setAbsMode(False)

def fitData(func,x,y,err,fields,fieldValues,idNums, nPar):
    global serialMode
    func.setXY(x,y)
    func.setIds(idNums)
    func.setErr(err)
    func.setFieldValues(fieldValues)
    func.setFields(fields)
    #print x
    #print y
    #print fieldValues
    #print idNums 
    func.setAbsMode(True)
    guesses = func.guess()
    boundaries = func.boundaries()
    sigma = []
    for bLow,bUp in zip(boundaries[0],boundaries[1]):
       sigma.append((bUp-bLow)/10.0)
    
    result = func.refine(guesses, boundaries[0], boundaries[1], sigma)
    fitPars = result.getPoint()

    if serialMode:
        simResult = func.simBounds(fitPars, boundaries[0], boundaries[1], sigma)
    else:
        simResult = func.simBoundsStream(fitPars, boundaries[0], boundaries[1], sigma)
    rex = func.getRex(fitPars)
    rexErr = func.getRexError()
    kex = func.getKex(fitPars)
    aicc = func.getAICc(fitPars)
    rms = func.getRMS(fitPars)

    return func.getParNames(),fitPars,simResult,aicc,rms,rex,rexErr,kex

def loadText(fileName):
    f1 = open(fileName,'rU')
    line = f1.readline()
    header =  line.strip().split('\t')
    useX = []
    peakRefs = []
    peakRef = 0
    for value in header[1:-2:2]:
        tau = float(value.strip())
        vcpmg = 1000.0/(2.0*tau)
        useX.append(vcpmg)
        peakRefs.append(str(peakRef))
        peakRef += 1

    resDataValues = {}
    peakNum = 0
    for line in f1:
        useY = []
        useErr = []
        values =  line.strip().split('\t')
        for (y,sdev) in zip(values[1:-2:2],values[2:-1:2]):
            useY.append(float(y.strip()))
            useErr.append(float(sdev.strip()))
        resNum = values[0].strip()
        resData = ResData(useX,useY,useErr,resNum,str(peakNum))
        resDataValues[int(resNum)] = resData
        peakNum += 1
    f1.close()
    return (resDataValues,peakRefs)



def loadMetaFile(fileName):
    f1 = open(fileName)
    for line in f1:
        if len(line) == 0:
            continue
        if line[0] == "#":
            continue
        values =line.split()
        if len(values) < 2:
            continue
        if values[0] == "field":
            fields = [float(values[1])]
        elif values[0] == "time":
            tCPMG = float(values[1])
        elif values[0] == "vcpmg":
            vcpmgs = [float(v) for v in values[2:]]
        else:
            print 'skip',line
    
    f1.close()
    return (tCPMG,vcpmgs,fields)

class CPMGExperiment:
    def __init__(self, field, temperature):
        self.field = field
        self.temperature = temperature

    def setData(self, resData):
        self.resData = resData

    def setPeakRefs(self, peakRefs):
        self.peakRefs = peakRefs

    def setTau(self, tau):
        self.tCPMG = tau

    def setV(self, vCPMG):
        self.vCPMG = vCPMG

class ResData:
    def __init__(self, x, y,err,  resNum,peakNum):
        self.x = x    
        self.y = y    
        self.err = err    
        self.resNum = resNum
        self.peakNum = peakNum

def loadDataFile(fileName):
    useXs = []
    useYs = []
    useErrs = []
    useFields = []
    resNums = []
    peakNums = []
    dupSDevs = []
    resDataValues = {}
    f1 = open(fileName)
    headers = f1.readline().strip().split('\t')
    firstIntensity = 3
    # skip first peakRef as it is for the reference which we don't use
    # maybe we should store it though
    peakRefs = [header[1:] for header in headers[firstIntensity+1:]]
    for line in f1:
        line = line.strip()
        #print line
        if line[0] == "#":
            continue
        values = line.strip().split('\t')
        if len(values) < 2:
            continue
        #print values
    #Peak	Residue	N	T1	T2	T11	T3	T4	T9	T5	T10	T12	T6	T7	T8
        peakNum = values[0].strip()
        resNum = values[1].strip()
        ref = float(values[firstIntensity].strip())
        useFieldValues = []
        useErr = []
        useX = []
        useY = []
        #print 'ref',ref
        if ref <= 0.0:
            continue
        #print ref,values[firstIntensity+1:]
        usedValues = [-1]*len(values[firstIntensity+1:])
        # fixme
        if float(values[firstIntensity+1]) > ref:
            print '#skip',values
            continue
        ok = True
        for i,v in enumerate(values[firstIntensity+1:]):
            if v != "NA":
                v = float(v)
                if v > 0.0:
                    if (v > ref):
                        ok = False
                        break
                    usedValues[i] = len(useX)
                    r2Eff = -math.log(v/ref)/tCPMG
                    useX.append(vcpmgs[i])
                    useY.append(r2Eff)
                    useFieldValues.append(fieldValues[i])
        if not ok:
            continue
        dupSDev = getSDevFromPeakDups(duplicates,usedValues,useY)
        useXs.append(useX)
        useYs.append(useY)
        useFields.append(useFieldValues)
        useErr = [dupSDev] * len(useX)
        resData = ResData(useX,useY,useErr,resNum,peakNum)
        useErrs.append(useErr)
        resNums.append(resNum)
        peakNums.append(peakNum)
        dupSDevs.append(dupSDev)
        resDataValues[int(resNum)] = resData
    
    f1.close()
    #return (resNums,useXs, useYs, useErrs, useFields, dupSDevs,peakNums,peakRefs)
    return (resDataValues,peakRefs)

class FitResult:
    def __init__(self, parNames, nGroupPars, nInGroup, fitPars, sdevs, aicc, rms, rex,rexErr, kex):
        self.parNames=parNames
        self.nGroupPars=nGroupPars
        self.aicc = aicc
        self.rms = rms
        self.rex = rex
        self.rexErr = rexErr
        self.kex = kex  
        self.values = {}
        self.sdevs = {}

        nNonGroup = len(parNames) - nGroupPars


        for i in range(nInGroup):
            for j in range(nGroupPars):
                parName = parNames[j]
                self.values[i,parName] = fitPars[j]
                self.sdevs[i,parName] = sdevs[j]

        for i in range(nInGroup):
            for j in range(nNonGroup):
                parName = parNames[nGroupPars+j]
                k = nGroupPars + i * nNonGroup + j
                self.values[i,parName] = fitPars[k]
                self.sdevs[i,parName] = sdevs[k]
        for i in range(nInGroup):
            if not "Rex" in parNames:
                self.values[i,"Rex"] = rex[i]
                self.sdevs[i,"Rex"] = rexErr[i]
        
def fitAllDataValuesOld(resDataValues,fields,peakRefs):
    for rD in resDataValues:
        aicMin = 1.0e9
        eqnLine = ""
        for eqn in ("NOEX","CPMGFAST","CPMGSLOW"):
            calcR.setEquation(eqn)
            parNames = calcR.getParNames()
            nPar = len(parNames)
            #print eqn,nPar
            parNames,fitPars,sdevs,aicc,rms,rex,kex =  fitData(calcR, rD.x, rD.y, rD.err, fields, rD.fields,nPar)
            eqnLine += "residue\t%s\tequation\t%8s\tRMS\t%6.2f\tAICc\t%6.2f" % (rD.resNum,eqn,rms,aicc)
            if aicc < aicMin:
                aicMin = aicc
                bestEqn = eqn
                bestRex = rex
                bestKex = kex
                bestRMS = rms
                bestR = fitPars[parNames.index('R2')]
                nPars = len(fitPars)
            eqnLine += "\tRex\t%7.2f\tKex\t%7.2f" % (rex,kex)
            for parName,fitPar,sdev in zip(parNames,fitPars,sdevs):
                eqnLine += "\t%3s\t%7.2f\t%s\t%7.2f" % (parName,fitPar,parName+'.err',sdev)
            eqnLine += '\n'
        resNum = rD.resNum
        outLine = "residue\t%s\tbest    \t%8s\tRMS\t%6.2f\tFRMS\t%6.2f\tR2\t%7.2f\tRex\t%7.2f\tKex\t%7.2f\n" % (resNum,bestEqn,rD.err[0],bestRMS,bestR,bestRex,bestKex)
        outLine += eqnLine
        outLine += "residue\t%s\txValues" % (resNum)
        for x in rD.x:
            outLine += "\t%6.2f" % (x)
        outLine += "\nresidue\t%s\tyValues" % (resNum)
        for y in rD.y:
            outLine += "\t%6.2f" % (y)
        outLine += "\nresidue\t%s\tpeakRefs" % (resNum)
        for peakRef in peakRefs:
            outLine += "\t%s" % (rD.peakNum+':'+peakRef)
        outLine += '\n'
        print outLine
        #alcR.dump(fitPars)
        #calcR.dump([16.6,22.5,1.8])

def writeHeader():
    eqnLine = "%s\t%5s\t%5s\t%5s\t%8s\t%6s\t%6s\t%s" % ("Residue","Peak","GrpSz","Group","Equation","RMS","AIC","Best")
    for parName in "R2","Rex","Kex","pA","dW":
        eqnLine += "\t%7s\t%7s" % (parName,parName+".sd")
    print eqnLine

def fitGroup(experiments,groupNums,groupID):

    x = []
    y = []
    f = []
    err = []
    idNums = []
    id = 0
    nInGroup = len(groupNums)
    for resNum in groupNums:
        fields = []
        for experiment in experiments:
            resDataValues = experiment.resData
            peakRefs = experiment.peakRefs
            field = experiment.field
            fields.append(experiment.field)
            rD = resDataValues[resNum]
            x += rD.x
            y += rD.y
            err += rD.err
            f += [field]*len(rD.x)
            idNums += [id]*len(rD.x)
        id += 1

    eqnLine = ""
    fitResults = {}
    aicMin = 1.0e9
       
    for eqn in ("NOEX","CPMGFAST","CPMGSLOW"):
        calcR.setEquation(eqn)
        parNames = calcR.getParNames()
        nGroupPars = calcR.getNGroupPars()
        nPar = len(parNames)
        parNames,fitPars,sdevs,aicc,rms,rex,rexErr, kex =  fitData(calcR, x, y, err, fields, f, idNums, nPar)
        fitResult  =  FitResult(parNames,nGroupPars, nInGroup, fitPars,sdevs,aicc,rms,rex,rexErr, kex)
        fitResults[eqn] = fitResult
        if fitResult.aicc < aicMin:
            bestEqn = eqn
            aicMin = aicc

    for id,resNum in enumerate(groupNums):
        rD = resDataValues[resNum]
        for eqn in ("NOEX","CPMGFAST","CPMGSLOW"):
            fR = fitResults[eqn]
            if eqn == bestEqn:
                best = "best"
            else:
                best = ""
            eqnLine = "%s\t%6s\t%5d\t%5d\t%8s\t%6.2f\t%6.2f\t%s" % (rD.resNum,rD.peakNum,nInGroup,groupID,eqn,fR.rms,fR.aicc,best)
            for parName in "R2","Rex","Kex","pA","dW":
                if (id,parName) in fR.values:
                    eqnLine += "\t%7.3f\t%7.3f" % (fR.values[id,parName],fR.sdevs[id,parName])
                else:
                    eqnLine += "\t\t"

            print eqnLine

def loadDataFiles(fileNames):
    global expList
    expList = []
    for fileName in fileNames:
        pathParts = os.path.split(fileName)
        fileTail = pathParts[-1][0:-4]

        parts = fileTail.split('_')
        temperature = float(parts[1])+273
        fieldValue = float(parts[2])
   
        (resDataValues,peakRefs) =  loadText(fileName)
        cpmgExp = CPMGExperiment(fieldValue, temperature)
        cpmgExp.setData(resDataValues)
        cpmgExp.setPeakRefs(peakRefs)
        expList.append(cpmgExp)

def loadProject(fileName):
    resProp = DataIO.loadParameters(fileName)
    return resProp

def getResidues(resProp):
    expDataSets = resProp.getExperimentData()
    resSet = set() 
    for expDataSet in expDataSets:
        residues = expDataSet.getResidues()
        for residue in residues:
           resSet.add(residue)
    residues = []
    for residue in resSet:
        residues.append([residue])
    return residues

def fitProject(resProp, groups, equationName):
    expDataSets = resProp.getExperimentData()
    residueFitter = ResidueFitter()
    if len(groups) == 0:
        groups = getResidues(resProp)
    groupID = 0
    for group in groups:
        print 'fit',group
        sgroup = [str(groupNum) for groupNum in group]
        resInfoList = residueFitter.fitResidues(resProp, sgroup, groupID, equationName)
        groupID += 1
        for resInfo in resInfoList:
            fitResNum = resInfo.resNum;
            resProp.addResidueInfo(str(fitResNum), resInfo)
    DataIO.saveParametersToFile('output.txt',resProp)

def fitGroups(groups):
    writeHeader()
    groupID = 0
    experiment = expList[0]
    groupDone = False
    for resNum in experiment.resData:
        inGroup = False
        for group in groups:
            if resNum in group:
                inGroup = True
                break
        if inGroup:
            if resNum == group[0]:
                fitGroup(expList,group,groupID)
                groupID += 1
        else:
            if not onlyGroups:
                fitGroup(expList,[resNum],groupID)
                groupID += 1

def parseArgs():
    global serialMode
    global onlyGroups
    serailMode = False
    parser = argparse.ArgumentParser(description="Analyze CPMG data")
    parser.add_argument("-s",dest='serialMode',action='store_true',help="Run calcutions serially (not using multiple cores")
    parser.add_argument("-e",dest='equationName',default=None,help="Fit specified equation, default is to fit all.")
    parser.add_argument("-o",dest='onlyGroups',action='store_true',help="Only fit data in groups")
    parser.add_argument("-g", dest="groupList",default='', help="Residues to fit in groups")
    parser.add_argument("-p", dest="projectFile",default='', help="Project file (.yaml) to load")
    parser.add_argument("fileNames",nargs="*")
    args = parser.parse_args()
    print args

    serialMode = args.serialMode
    onlyGroups = args.onlyGroups
    fileNames = args.fileNames
    projectFile = args.projectFile
    equationName = args.equationName
    print 'fileNames',fileNames

    groups = [] 
    print 'gL',args.groupList
    if len(args.groupList) > 0:
        sgroups = args.groupList.split(' ')
        print 'sg',sgroups
        groups = []
        for group in sgroups:
            group = [int(groupNum) for groupNum in group.split(',')]
            groups.append(group)
        print 'g',groups

    if projectFile != '':
        resProp = loadProject(projectFile)
        fitProject(resProp, groups, equationName)
    else:
        loadDataFiles(fileNames)
        fitGroups(groups)

if len(sys.argv) > 1:
    arg0 = sys.argv[1]
    if arg0.endswith(".py"):
        sys.argv.pop(0)
        runpy.run_path(arg0)
    else:
        sys.argv[0] = "cpmgcalc.py" 
        parseArgs()

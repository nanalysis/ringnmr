def generateNoise(mean, std, clean_ann_input):
    """
    mean (float|int)
    std (float)
    clean_ann_input (list<float>)
    """

    import random as rd
    return [(value+rd.gauss(mean, std)) for value in clean_ann_input]

def getSavedNeuralNetwork(annDir):
    import os
    # FIXME : path should not be from personal comp
    savedANNPath = "/Users/teddycolon/Desktop/ASRC/ojAlgo-ANN/{}/saved-networks/".format(annDir)
    savedNetworks = os.listdir(savedANNPath)

    if savedNetworks:
        savedNetworkFile = ""
        for elem in savedNetworks:
            if elem.endswith(".txt"):
                savedNetworkFile = elem
                break
        if savedNetworkFile:
            return (savedANNPath + savedNetworkFile)
        else:
            raise Exception("Could not find a saved ANN file.")
    else:
        raise Exception("The directory '../saved-networks' is empty.")

def generateOutputString(annInputInfoList, targetPars):
    """
    annInputInfoList (list<floats>) : ANN input values obtained from 'organizeANNInfo(..)'
    targetPars (list<floats>) : randomized target parameters used to calculate Y-Axis intensity values
                                for the CEST r1RhoPerturbation profile

    returns _ (string) : string containing ANN input and target value information 
    """
    inputStringBuilder = [str(val) for val in annInputInfoList]
    annTargetPars = [str(val) for val in targetPars]
    annInputString = ','.join(inputStringBuilder)
    annTargetParsString = ','.join(annTargetPars)
    return ':'.join([annInputString, annTargetParsString])


def getSimX(nPoints, lowBound, uppBound):
    """
    nPoints (int) : number of points describing CEST profile
    lowBound (float) : lower offset bound (x axis)
    uppBound (float) : upper offset bound (x axis)

    returns xVals (list) : list of x values in range of low/upper bounds
    """
    xVals = []
    delta = (uppBound - lowBound) / (nPoints + 1)
    value = lowBound
    for _ in range(nPoints):
        xVals.append(value)
        value += delta
    return xVals

def getMean(values):
    """
    values (list<floats|ints>) : list of floats to find average

    return _ (float) : the average of the values in the list
    """
    if isinstance(values, list):
        if len(values) > 1:
            if isinstance(values[0], float) or isinstance(values[0], int):
                sumValues = 0.0
                for elem in values:
                    sumValues += elem
                return (sumValues / len(values))
            elif isinstance(values[0], list) or isinstance(values[0], tuple):
                firstRow = values[0]
                allSameSize = True
                for i, _ in enumerate(values):
                    nextI = i + 1
                    if nextI < len(values):
                        if len(firstRow) != len(values[nextI]):
                            allSameSize = False
                if allSameSize:
                    zipResults = list(zip(*values))
                    finalResult = []
                    for result in zipResults:
                        finalResult.append(getMean(list(result)))
                    return finalResult
                else:
                    raise ValueError("Each sublist must be the same size.")
        elif len(values) == 1:
            return values
        else: 
            raise ValueError("The input list is empty!")
    else:
        raise TypeError("The type of input '{}' is not valid!".format(values))


def scaleValue(value, maxValue, minValue):
    """
    value (float) : the value to scale
    maxValue (float) : maximum value
    minValue (float) : minimum value

    return _ (float) : value b/t 0 - 1
    """
    return (value - minValue) / (maxValue - minValue)

def revertScale(value, maxValue, minValue):
    """
    value (float) : value scaled from 0 - 1
    maxValue (float) : maximum value used to scale the scaled value
    minValue (float) : minimum value used to scale the scaled value

    return _ (float) : original value 
    """
    return (((value*maxValue) - (value*minValue)) + minValue)

def closeFile(fileObj):
    """
    fileObj (file object) : An open file object

    returns _ (file object): closed file object
    """
    return fileObj.close()

def makeFile(fileName):
    """
    fileName (String) : The name of the file to store train/test examples

    returns _ (file object): file object with the name provided
    """
    return open(fileName, 'w')


def getSimX(nPoints, lowBound, uppBound):
    """
    nPoints (int) : number of points describing CEST profile
    lowBound (float) : lower offset bound (x axis)
    uppBound (float) : upper offset bound (x axis)

    returns xVals (list) : list of x values in range of low/upper bounds
    """
    xVals = []
    delta = (uppBound - lowBound) / (nPoints + 1)
    value = lowBound
    for _ in range(nPoints):
        xVals.append(value)
        value += delta
    return xVals


def trainNetwork(trainFileName, validationFileName, nInputNodes, nOutputNodes, saveNetworkDir, scaleValuesString):
    """
    Purpose :
        retrieve training information from file provided and train/save the neural network. 

    fileName (string) : name of file that contains training data

    return _ : saves network into the directory specified.
    """
    import os
    from ANN import CESTTrain as ct
    from ANN import NeuralNetworkUtils as netutils
    
    if os.path.isfile(validationFileName):
        CESTtraining = ct(trainFileName, validationFileName)
    else:
        CESTtraining = ct(trainFileName)
    
    CESTtraining.setScales(scaleValuesString) # setting scale values to save after training the neural network

    result = CESTtraining.cestR1RhoPANNRun("R1RHOPERTURBATION", nInputNodes, nOutputNodes)

    if result:
        trainedNetwork = result.getFirst()
        networkInfo = result.getSecond()

        savedNeuralNetworkFile = netutils.saveNeuralNetwork(saveNetworkDir, trainedNetwork, networkInfo)
    else:
        raise Exception("Result from training neural network run is None/Null.")

def updateProgress(progress, total):
    """
    progress (int) : current value
    total (int) : total values
    """
    import time
    import sys
    time.sleep(0.005)
    sys.stdout.write("\r %d out of %d" % (progress, total))
    sys.stdout.flush()


def getRMSE(actual, expected):
    from math import sqrt
    sumSq = 0.0
    nPars = len(actual) if len(actual) == len(expected) else -1
    if nPars < 0:
        raise AssertionError("Size {0} != {1}".format(len(actual), len(expected)))
    diff = 0.0
    for i,_ in enumerate(actual):
        diff = actual[i] - expected[i]
        sumSq += diff ** 2

    return sqrt(sumSq/nPars)

def displayData(*dataInfo):
    if len(dataInfo) == 5:
        X, Y, plotTitle, plotXLabel, plotYLabel = dataInfo
        X = [float(xStringVal) for xStringVal in X.split(',')]
        Y = [float(yStringVal) for yStringVal in Y.split(',')]
        
        import matplotlib.pyplot as plt

        plt.scatter(X,Y,color='r')
        plt.title(plotTitle)
        plt.xlabel(plotXLabel)
        plt.ylabel(plotYLabel)
        plt.show()

def displayCestProfiles(xValues, yValues, counter):
    import matplotlib.pyplot as plt
    
    x = [float(elem) for elem in xValues.split(',')]
    y = [float(elem) for elem in yValues.split(',')]
    xLower, xUpper = min(x) - 1, max(x) + 1


    ax.scatter(x,y,marker='.')
    ax.plot(x,y,marker='o')
    ax.set(xlim=[xLower, xUpper], yLim=[0.0, 1.0])
    ax.set(title="CEST Profile Simulated Data {}".format(count), xlabel="offset (ppm)", ylabel="Intensity (I/I0)")

    plt.show()

def comparisonFile(labelList, actual, iterative, network, compareFile):
    if len(actual) == len(iterative) == len(network) == len(labelList):
        #labelList = ["KEX", "pB", "deltaA", "deltaB", "r1A", "r1B", "r2A", "r2B"]
        #labelList = ["KEX", "pB", "r2A", "r2B"]
        refString = "\t%-10s <---> %-10s <---> %-10s\n" % ("actual", "current", "network")
        compareFile.write(refString)
        for i,_ in enumerate(actual):
            string = "%-6s: " % (labelList[i])
            string += "%-10f <---> %-10f <---> %-10f\n" % (actual[i], iterative[i], network[i])
            compareFile.write(string)
        compareFile.write("+++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
    else:
        closeFile(compareFile)
        print("actual : {0},\niterative : {1},\nnetwork : {2}".format(actual, iterative,network))
        raise AssertionError("The length of actual, iterative, network arrays are not equal. ({0} != {1} != {2})".format(len(actual), len(iterative), len(network)))

def displayRMSEBarPlot(**kwargs):

    experimentType = kwargs.get("expType")
    mode = kwargs.get("mode")
    nTests = kwargs.get("nTests")
    rmseListString = kwargs.get("rmseList")
    nFields = kwargs.get("nFields")
    plotTitle = kwargs.get("title", "")
    saveFileName = kwargs.get("saveFN", '')

    import matplotlib.pyplot as plt
    import numpy as np

    if nTests and rmseList:
        nTests = int(nTests)
        rmseList = [[float(val) for val in pair.split(':') if val] for pair in rmseListString.split(',') if pair]
    else:
        raise ValueError("Values for 'nTests' and/or 'rmseListString' do not exist.")

    currGuess = []
    annGuess = []
    for i in range(nTests):
        currRMSE, annRMSE = rmseList[i]
        currGuess.append(currRMSE)
        annGuess.append(annRMSE)
    
    ind = np.arange(nTests)
    width = 0.35
    plt.bar(ind, currGuess, width, label='Iterative Guess')
    plt.bar(ind + width, annGuess, width, label='Network Guess')

    plt.xlabel('Number of Tests')
    plt.ylabel('RMSE')

    if mode and experimentType.startswith("cpmg") and nFields:
        fieldListString = kwargs.get("fieldList")
        if fieldListString:
            fields = [float(val) for val in fieldListString.split(',') if val]
        else:
            raise ValueError("Value for 'fieldListString' does not exist.")
        plt.title('CPMG {0} for {1} field: RMSE Scores (FIELDS: {2}-{3})'.format(mode, nFields, min(fields), max(fields)))
    else:
        plt.title(plotTitle)

    plt.xticks(ind + width / 2, [str(i) for i in range(1,nTests+1)])
    plt.legend(loc='best')
    if saveFileName:
        plt.savefig(saveFileName)
    plt.show()


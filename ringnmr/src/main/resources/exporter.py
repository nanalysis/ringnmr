import random
import java
class Plot():
    # Represents an entire subplot, contains all the data necessary
    def __init__(self, xlabel, ylabel, title, type, rawData = None, fittedData = None):
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.title = title
        self.plotType = type # will be either barplot or scatter and will effect plot of rawData
        # To store dataset
        self.x = []
        self.y = []
        self.error = []
        self.fittedX = []
        self.fittedY = []

        if rawData:
            self.unpackData(rawData,fittedData)

        # To store the variable names of the data
        self.xVar = []
        self.yVar = []
        self.errorVar = []
        self.fittedXVar = []
        self.fittedYVar = []

    def __str__(self):
        lines = '/n'.join(['xlabel: ' + self.xlabel, 'ylabel: ' + self.ylabel, 'title: ' + self.title, 'plotType: ' + self.plotType])
        return lines

    def __getitem__(self,item):
        returns = {'xVar': self.xVar, 'yVar': self.yVar, 'errorVar': self.errorVar, 'fittedXVar': self.fittedXVar, 'fittedYVar':self.fittedYVar, 'x': self.x, 'y': self.y, 'error': self.error, 'fittedX': self.fittedX, 'fittedY': self.fittedY}
        return returns[item]

    def setColor(self, color):
        self.color = color # can be an index for a variable or char

    def unpackData(self, rawData, fittedData=None, varDict = None):
        x,y, error = table(rawData)
        self.x.append(x)
        self.y.append(y)
        self.error.append(error)

        if fittedData:
            x,y = table(fittedData)
            self.fittedX.append(x)
            self.fittedY.append(y)

        if varDict: # varDict to contain info about var, baseName and optional index
            names = ['X','Y','Er']
            baseName = varDict['baseName']
            rawVarNames = []
            for name in names:
                varNameSequence = [baseName, 'raw',name,varDict['index']] if 'index' in varDict else [baseName,'raw',name]
                rawVarNames.append('_'.join(varNameSequence))
            fittedVarNames = []
            if fittedData:
                for name in names[:-1]:
                    varNameSequence = [baseName,'fitted',name,varDict['index']] if 'index' in varDict else [baseName,'fitted',name]
                    fittedVarNames.append('_'.join(varNameSequence))
            self.addDataVars(rawVarNames, fittedVarNames)

    def addDataVars(self, rawVarNames, fittedVarNames=None):
        # Called to store variable names for data
        self.xVar.append(rawVarNames[0])
        self.yVar.append(rawVarNames[1])
        self.errorVar.append(rawVarNames[2])
        if fittedVarNames:
            self.fittedXVar.append(fittedVarNames[0])
            self.fittedYVar.append(fittedVarNames[1])

    def setLegend(self, legend):
        self.legend = legend # should be a list of strings in order of the variables plotted

    def setRanges(self, xRange, yRange):
        self.xrange = xRange
        self.yrange = yRange


class Writer():
    def __init__(self, configData, data, barChartData):
        # Initializes and pythonizes the data for the Writer
        self.configData = configData;
        # The config data is primarily for the residue charts

        self.figures = [] # to contain a list of Plot objects

        #data and barChartData is all data for the two figures to be generated
        self.data = data;
        self.barChartData = barChartData;

        # Turns the java objects into lists and dicts
        self.pythonizeExportData()
        self.type = self.configData['exportType']

        #if true will scale color by 255
        self.scaleColor = (self.type == 'grace')

        #if true will create variables for the data
        self.doVariables = (self.type != 'grace') # Grace does not have variables

        # Creates list of colors for the object
        self.setRGBColors()

        self.template = getFormatData(self.type);
        self.lines = [] # Contain the lines that will get printed out to page
        # Vars will include title, xlabel, ylabel, xlims, ylims, dataSet

    def pythonizeExportData(self):
        # TODO fixup pythonize.py
        self.configData = dict(self.configData)
        self.configData['ranges'] = list(self.configData['ranges'])
        self.data = pythonizeChart(self.data)
        self.barChartData = pythonizeChart(self.barChartData)

    def writeFromExportData(self):
        self.addTop()       #Sets up file header but also adds color variables
        self.configureLayout() #Generates residueIndices to point to which data should be on each subplot and finds nRows, nCols
        if self.doVariables:
            self.setUpPlots()   #Adds list of plot objects to figures
        else:
            self.addLabels()
        self.writeData()
        if 'endplot' in self.template:
            self.lines.append(self.template['endplot'])
        self.writeToFile()

    def setUpPlots(self):
        figure = []
        baseTitle = self.configData['title']
        xlabel = self.configData['xlabel']
        ylabel = self.configData['ylabel']
        ranges = self.configData['ranges']
        for res in self.residueIndices:
            title = ' '.join([baseTitle,'for res',res])
            subplot = Plot(xlabel, ylabel, title, 'scatter')
            legend = []
            colors = []
            for i, index in enumerate(self.residueIndices[res]):
                colors.append(min(index+1,len(self.colors))+1)
                plotData = self.data[index]
                graphTitle = plotData['graphTitle']
                rawData = plotData['rawData']
                fittedData = plotData['fittedData']
                legend.append(graphTitle.split(':')[0]) # this has all the information of the plot
                varDict = {'baseName' : 'res'+str(res), 'index' : str(i+1)}
                subplot.unpackData(rawData,fittedData,varDict)
            subplot.setColor(colors)
            subplot.setLegend(legend)
            subplot.setRanges(ranges[:2],ranges[2:])
            figure.append(subplot)
        self.figures.append(figure)
        figure = []
        for plot in self.barChartData:
            subplots = []
            title = ''
            xlabel = 'Residue Number'
            ylabel = ''
            ranges = plot['ranges']
            del plot['ranges']
            legend = []
            subplot = None
            colors = ['green','blue','red'] if self.type == 'r' else ['g','b','r']
            colors = colors[0:min(len(plot),3)] # Colors of chart
            for i, plotData in enumerate(plot):  #plot is a dictionary with keys of plot names and values of data values
                genTitle, fit, _, dataFitted = plotData.split("|") # these keys have title, type of fit, _ and type of data fitted
                if not title:
                    title = genTitle + ":" + dataFitted
                    ylabel = dataFitted
                    subplot = Plot(xlabel, ylabel, title, 'barplot')
                legend.append(fit)
                varDict = {'baseName': '_'.join([dataFitted,fit])}
                subplot.unpackData(plot[plotData],varDict=varDict)
            subplot.setColor(colors)
            subplot.setRanges(ranges[:2],ranges[2:])
            subplot.setLegend(legend)
            figure.append(subplot)
        self.figures.append(figure)

    def writeFigures(self):
        subplotTemplate = self.template['subplot']
        for index, figure in enumerate(self.figures):
            if index: #index == 1
                self.nRows = len(figure); self.nCols = 1;
                self.lines.append(self.template['newplot'])
            plotVars = []
            for i,subplot in enumerate(figure):
                if self.type == 'python':
                    line = subplotTemplate.format(self.nRows,self.nCols, i+1)
                    self.lines.append(line)
                plot = self.writeSubplot(subplot, i, returnStr=self.type=="r",)
                # Below will only occur if writeSubplot returns something
                if plot:
                    template = self.template['var']
                    varName = '_'.join(['plot',str(i)])
                    line = template.format(plot, varName)
                    self.lines.append(line)
                    plotVars.append(varName)
            if self.type == 'r':
                line = subplotTemplate.format(','.join(plotVars), self.nRows, self.nCols)
                self.lines.append(line)


    def writeSubplot(self, subplot, index, returnStr=False):
        template = self.template['var']
        datasetNames = ('x','y','error','fittedX','fittedY') if subplot.fittedX else ('x','y','error')
        subplotData = []
        subplotVars = []
        for datasetName in datasetNames:
            subplotData.append(subplot[datasetName])
            subplotVars.append(subplot[datasetName+'Var'])
        for i, _ in enumerate(subplotData):
            for j, _ in enumerate(subplotData[i]):
                varName = subplotVars[i][j]
                varData = subplotData[i][j]
                varData = self.createArr(varData)
                line = template.format(varData,varName)
                self.lines.append(line)
        plotComponents = [] # Instead of making individual lines, just find components. This is mostly due to R (subplots must have + between them)
        if self.type == "r":
            dataFrames = self.createDataFrames(subplot, index)
            rawDF = dataFrames['rawDF']  # variable for scatter/barcharts
            fittedDF = dataFrames['fittedDF'] # variable for line graphs
            plotComponents.append(self.template['firstplt'].format(rawDF))
            plotComponents.append(self.template[subplot.plotType])
            if fittedDF:
                plotComponents.append(self.template['line'].format(fittedDF))
            plotComponents.append(self.template['error'])
        else:
            for i, _ in enumerate(subplotVars[0]):
                color = subplot.color[i] if subplot.plotType == 'scatter' else subplot.color[i%3]
                xyEr = [subplotVars[0][i], subplotVars[1][i], subplotVars[2][i], color]
                if subplot.plotType == 'barplot':
                    numBars = len(subplotVars[0])
                    if numBars%2==0:
                        shifts = [j+0.5 for j in range(-(numBars/2),0)] + [j-0.5 for j in range(1,numBars/2+1)]
                    else:
                        shifts = range(-(numBars/2),numBars/2+1)
                    xyEr.append(shifts[i])
                xyEr.append(subplot.legend[i])
                plotComponents.append(self.template[subplot.plotType].format(*xyEr))
                if len(subplotVars) > 3:
                    fittedXY = [subplotVars[3][i],subplotVars[4][i], color]
                    plotComponents.append(self.template['line'].format(*fittedXY))

        configComponents = self.createConfigLines(subplot, numSubPlot = index, numPlots = len(subplotVars[0]))
        components = plotComponents + configComponents
        line = ' + '.join(components) if self.type == 'r' else '\n'.join(components)
        if returnStr:
            return line
        else:
            self.lines.append(line)

    def createDataFrames(self, subplot, index):
        dataFrameVar = self.createDataFrame(subplot,('x','y','error'),False, index)
        fittedDataFrameVar = ''
        if subplot.fittedX:
            fittedDataFrameVar = self.createDataFrame(subplot,('x','y'), True, index)
        return {'rawDF': dataFrameVar, 'fittedDF': fittedDataFrameVar}

    def createDataFrame(self, subplot, varSets, fitted, index):
        templateVar=self.template['var']
        for varSet in varSets:
            if fitted:
                line = templateVar.format(self.createArr(subplot['fitted'+varSet.upper()+'Var']),varSet)
            else:
                line = templateVar.format(self.createArr(subplot[varSet+'Var']),varSet)
            self.lines.append(line)
        testData = subplot.fittedX if fitted else subplot.x
        groups = self.createArr(createGroupings(testData,subplot.legend),True)
        self.lines.append(templateVar.format(groups,'group'))
        template = self.template['frame']
        frameArr = varSets + tuple(['group'])
        frameInput = ','.join(frameArr)
        frame = template.format(frameInput)
        dataFrameVar = '_'.join(['df', 'fitted' if fitted else 'raw', str(index)])
        self.lines.append(templateVar.format(frame,dataFrameVar))
        return dataFrameVar


    def createConfigLines(self,subplot, numSubPlot = None, numPlots = None):
        templates = ['title','xlabel','ylabel','xlim','ylim']
        values = [subplot.title, subplot.xlabel, subplot.ylabel, self.createArr(subplot.xrange), self.createArr(subplot.yrange)]
        components = [self.template[templates[i]].format(values[i]) for i,val in enumerate(templates)]
        if self.type == 'r':
            template = self.template['pltcolor']
            colors = subplot.color
            if numPlots > 3 and subplot.plotType == 'barplot': # 3 different colors
                for i in range(numPlots-3):
                    colors.append(colors[i%3])
            if subplot.plotType == 'scatter':
                colors = ['color'+str(color) for color in colors]
            colors = self.createArr(colors) if subplot.plotType =='scatter' else self.createArr(colors, isString=True)
            components.append(template.format(colors))
        else:
            if subplot.plotType == 'barplot':
                template = self.template['xtick']
                components.append(template.format(int(subplot.xrange[0]),int(subplot.xrange[1])))
            components.append(self.template['legend'])
        return components

    def configureLayout(self):
        import math
        residues = [plot['graphTitle'].split(':')[1] for plot in self.data] #makes list of
        residueIndices = {}
        for i, val in enumerate(residues):
            residueIndices.setdefault(val, [])
            residueIndices[val].append(i)
        self.residueIndices = residueIndices
        numPlots = len(residueIndices)
        for i in range(9):        # Caps the num of subplots at 81
            if i*i >= numPlots:
                nCols = i;
                break
        nRows = int(math.ceil(float(numPlots)/nCols))
        self.nCols = nCols
        self.nRows = nRows

    def addTop(self):
        if 'encode' in self.template:
            self.lines.append(self.template['encode'])
        if 'import' in self.template:
            self.lines.append(self.template['import'])

        if 'colors' in self.configData:
            for i in enumerate(self.colors):
                # 0 and 1 for background colors in grace
                rgb = self.colors[i[0]]
                line = self.template['color'].format(i[0]+2,rgb[0], rgb[1], rgb[2])
                self.lines.append(line)

    def createArr(self, arr, isString=False):
        template = self.template['arr']
        arr = [str(x) for x in arr]
        if isString:
            arr = ['"'+datum+'"' for datum in arr]
        arrStr = ','.join(arr)
        line = template.format(arrStr)
        return line

    def writeData(self):
        if self.doVariables:
            self.writeFigures()
        else:
            self.prepareGraceData()
            self.writeGraceData()

    def prepareGraceData(self):
        datasets = self.data
        keys = ['rawData','fittedData']
        graphNum = 0;
        for dataset in datasets:
            for key in keys:
                plotType = 'scatter' if key == 'rawData' else 'line'
                lines = list(self.template[plotType])
                self.lines += [line.format(graphNum,graphNum/2 + 2) for line in lines]
                graphNum += 1

    def writeGraceData(self):
        datasets = self.data
        keys = ['rawData','fittedData']
        for dataset in datasets:
            for key in keys:
                type = 'xydy' if key == 'rawData' else 'xy'
                line = self.template['dataHeader'].format(type)
                self.lines.append(line)
                data = dataset[key]
                data = ['\t'.join([str(val) for val in data[i]]) for i,_ in enumerate(data)]
                self.lines += data
                self.lines.append('&')

    def setRGBColors(self):
        self.colors = [[color.getRed(), color.getGreen(), color.getBlue()] for color in self.configData['colors']]
        if self.scaleColor:
            self.colors = [[int(hue*255) for hue in color] for color in self.colors]

    def writeToFile(self):
        with open(self.configData['file'], 'w') as target:
            writeOut = '\n'.join(self.lines).encode('utf8')
            target.write(writeOut)

    def addLabels(self):
    # Only used for R
        iterableItems = ('title','xlabel','ylabel','ranges') # Ranges is xmin, xmax, ymin ymax
        for item in iterableItems:
            template = self.template[item]
            value = self.configData[item]
            if item == 'ranges':
                line = template.format(*value)
            else:
                value = resolveUnicode(value,self.type)
                line = template.format(value,item)
        self.lines.append(line)


def resolveUnicode(value, type):
    import unicodedata
    if type =='python':
        value = repr(value)
    elif type =='r':
        value = '"' + value + '"'
    else:
        # TODO setup a lookup table of sorts, to convert the unrecognizable Unicode into other values
        value = unicodedata.normalize('NFKD',value).encode('ascii','ignore')
    return value


def getFormatData(type):
    from org.yaml.snakeyaml import Yaml
    from java.io import ByteArrayInputStream
    formatData = loadResource('template.yaml')
    input = java.io.ByteArrayInputStream(formatData)
    yaml = Yaml()
    data = yaml.load(input)
    return data[type]


def loadResource(resourceName):
    from java.lang import ClassLoader
    from java.io import BufferedReader
    from java.io import InputStreamReader
    cl = ClassLoader.getSystemClassLoader()
    istream = cl.getResourceAsStream(resourceName)
    lines = ""
    if istream == None:
        raise Exception("Cannot find '" + resourceName + "' on classpath")
    else:
        reader = InputStreamReader(istream)
        breader = BufferedReader(reader)
        while True:
            line = breader.readLine()
            if line == None:
                break
            if lines != '':
                lines += '\n'
            lines += line
        breader.close()
    return lines


def table(data):
    data = [line.split() for line in data] if isinstance(data[0],basestring) else data
    tabled = []
    for i in range(len(data[0])):
       col = [x[i] for x in data]
       tabled.append(col)
    return tabled


def pythonizeChart(chart):
    chart = list(chart)
    for i, subplot in enumerate(chart):
        chart[i] = dict(subplot)
        for key in subplot:
            if key == 'graphTitle':
                continue
            chart[i][key] = list(subplot[key])
            if key == 'ranges':
                continue
            for j, point in enumerate(subplot[key]):
                chart[i][key][j] = list(point)
    return chart


def createGroupings(sampleData, legend):
    groups = []
    for i, label in enumerate(legend):
        groups += [label]*len(sampleData[i])
    return groups

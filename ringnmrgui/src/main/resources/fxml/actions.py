from java.lang import System
from java.util import ArrayList
from javafx.application import Platform
from javafx.stage import FileChooser
from javafx.stage.FileChooser import ExtensionFilter
from org.comdnmr.gui import ChartUtil
from org.comdnmr.gui import SecondaryStructure
from org.comdnmr.gui import PyController

def drawChart(event):
    node = event.getSource()
    chartNode = ChartUtil.findNode(node.getScene(),"sschart")
    if chartNode == None:
        print 'no chart node'
    else:
        xValues = [0.0,1.0,2.0,3.0,4.0]
        yValues = [2.0,1.0,0.5,3.0,2.0]
        data = ChartUtil.makeChartSeries(xValues,yValues)
        chartNode.setData(data)

def pushMe(event):
    print 'pyaction push from actions.py'
    print event
    node = event.getSource()
    print node
    print Platform.isFxApplicationThread()
    loadFile(event) 
    #System.exit(1)

def clearChart(event):
        node = event.getSource()
        chartNode = PyController.mainController.getActiveChart()
        if chartNode == None:
            print 'no chart node'
            return

        chartNode.getData().clear()

def printChart(event):
        node = event.getSource()
        chartNode = PyController.mainController.getActiveChart()
        if chartNode == None:
            print 'no chart node'
            return

        ChartUtil.print(chartNode)

def loadFile(event):
        node = event.getSource()
        chartNode = PyController.mainController.getActiveChart()
        if chartNode == None:
            print 'no chart node'
            return
        fileChooser = FileChooser()
        fileChooser.setTitle("Open FID")
        #filter = ExtensionFilter("All Files", ["*"])
        #fileChooser.getExtensionFilters().add(filter)
        print node
        selectedFile = fileChooser.showOpenDialog(node.getScene().getWindow())
        if selectedFile != None:
            print selectedFile
            data = ChartUtil.loadChartData(selectedFile.toString())
            print 'node',chartNode
            print 'data',data
            print 'cdata',chartNode.getData()
            chartNode.setUnifyYAxes(False)
            chartNode.getData().addAll(data)

def addChart(event):
    node = event.getSource()
    newChart = PyController.mainController.addChart()

def loadCPMGFile(event):
        node = event.getSource()
        fileChooser = FileChooser()
        fileChooser.setTitle("Open Chart Data")
        #filter = ExtensionFilter("All Files", ["*"])
        #fileChooser.getExtensionFilters().add(filter)
        print node
        selectedFile = fileChooser.showOpenDialog(node.getScene().getWindow())
        if selectedFile != None:
            nMaps = ChartUtil.getNMaps()
            nodeName = "sschart"
            ChartUtil.loadParameters(selectedFile.toString())



def loadSecondaryStructure(event):
    node = event.getSource()
    #chartNode = ChartUtil.findNode(node.getScene(),"xyBarChart0")
    chartNode = PyController.mainController.getActiveChart()
    ssRegion = ChartUtil.findNode(chartNode.getScene(),"ssregion")
    if chartNode == None:
        print 'no chart node'
        return
    fileChooser = FileChooser()
    fileChooser.setTitle("Open Secondary Structure")
    selectedFile = fileChooser.showOpenDialog(chartNode.getScene().getWindow())
    if selectedFile != None:
        f1 = open(selectedFile.toString(),'r')
        ssValues = ArrayList()
        for line in f1:
            line = line.strip()
            if len(line) == 0 or line[0] == '#':
                continue
            fields = line.strip().split()
            start = int(fields[0])
            end = int(fields[1])
            type = fields[2]
            label = ''
            if (len(fields) > 3):
                label = fields[3]
            color = 'gray';
            if (len(fields) > 4):
                color = fields[4]
            ss = SecondaryStructure(start, end, type, label, color)
            ssValues.add(ss)
        f1.close()
        #chartNode.addSS(ssValues)
        ssRegion.setChart(chartNode, ssValues)

def setSecondaryStructure(chartNode):
    ss = SecondaryStructure(10,20,0,'A')
    ssValues = ArrayList()
    ssValues.add(ss)
    #chartNode.addCanvas()
    chartNode.addSS(ssValues)

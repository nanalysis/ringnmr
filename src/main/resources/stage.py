from javafx.application import Platform
from javafx.stage import Stage
from javafx.scene import Scene
from javafx.fxml import FXMLLoader
from java.net import URL
from org.comdnmr.fit.gui import MainApp
from java.lang import ClassLoader
from sys import argv


class PyStage:
    stage = None
    scene = None
    controller = None
    def __init__(self, fxml=None, stage=None):
        if stage == None:
            self.stage = Stage()
        else:
            self.stage = stage
        if fxml != None:
            self.load(fxml)

    def load(self,urlString):

        url = ClassLoader.getSystemClassLoader().getResource(urlString);
        #url = URL(urlString)
        loader = FXMLLoader(url)
        root = loader.load()
        self.controller = loader.getController()
        self.scene = Scene(root)
        self.stage.setScene(self.scene)
        self.stage.show()

    def title(self,title):
        self.stage.setTitle(title)

def onActionOff(nodeID):
    actionFunctions[nodeID]()
    #print 'action on',str(nodeID)

def onAction(node):
    print "onaction"
    #actionFunctions[nodeID]()
    #print 'action on',str(nodeID)

def pushme():
    print 'pushed'

actionFunctions={}

actionFunctions['pushme'] = pushme

#pystage = PyStage("file:src/main/resources/fxml/NMRScene.fxml",MainApp.primaryStage)
pystage = PyStage("fxml/NMRScene.fxml",MainApp.primaryStage)
pystage.title("CoMD/NMR Dynamics Analysis")


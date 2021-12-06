import random

class Writer(): 

  def __init__(self):
    # The parameters to be included
    self.parameters = []

    self.subtitle = ""
    self.title = ""
    self.strings = [] # To contain strings to add to the view. Then to be added to params after all set to correct view
    self.top = .9; # Give room to the title. Will be value of y top position of graph.
    self.bottom = .1; # just provide some cushion. Will be value of y bottom position of graph. (each string will raise)
    self.stringCushion = .035; # spacing between the strings on view

    self.dataSets = [] # Each dataset is to be a dictionary: ['data'] to contain the raw data, ['plot'] to be scatter or line or (soon hopefully bar)

    self.axes = {'label' : ['"XLabel"','"YLabel"'], 'min' : [], 'max' : [], 'minorTick' : [], 'majorTick' : []};

    self.dataGroupings = []
    # The following should be sets of lists of different data. 
    self.dataGroups = []

  ### These Params are to be set before creating everything ###
  def setTitle(self, title):
    self.title = title
  
  def setSubtitle(self, subtitle):
    self.subtitle = subtitle
    self.top = .85

  def setString(self, string):
    string = '"' + string + '"'
    string = {'string': string, 'position': self.bottom}
    self.strings.append(string);
    self.bottom += self.stringCushion; 

  def setAxes(self, axes):
    for key in axes:
      self.axes[key] = axes[key]
    
  def addDataSets(self, dataSets):
    for dataSet in dataSets:
      self.dataSets.append(dataSet)
 
  def addDataSet(self, dataSet):
    self.dataSets.append(dataSet)

  def setStrings(self, strings):
    i = 0;
    for string in strings:
      self.setString(string)
      i += 1;
      if i == len(strings):
        self.bottom += .05
  
  def addRawData(self, data, color=-1):
    dataSet = {'data': data, 'plot': 'scatter'}
    if color >= 0:
      dataSet['color']=color
    self.addDataSet(dataSet)

  def addFittedData(self, data, color=-1):
    dataSet = {'data':data, 'plot':'line'}
    if color >= 0:
      dataSet['color']=color
    self.addDataSet(dataSet)

  ### These Params are to be set before creating everything ###

  #RECALL data set is a dictionary with 'data' and 'plot' where 'plot' is either scatter or line

  def addParam(self, param, value):
    line = " ".join(['@', param, value]) + '\n'
    self.parameters.append(line)

  def addBulkParams(self, params, values):
    for i in range(len(params)):
      self.addParam(params[i],values[i])

  def addTitleHead(self):
    isSubtitle = (self.subtitle != "") 
    isTitle = (self.title != "")
    
    if isTitle:
      title = '"' + self.title + '"'
    else:
      title = '"Exported Graph"'
      
    params = ['CLEAR','TITLE']
    values = ['STRING',title]
    if isSubtitle:
      subtile = '"' + self.subtitle +'"'
      params.append('SUBTITLE'); 
      values.append(self.subtitle);
    self.addBulkParams(params, values)

  def addViews(self):  # to be done after all the strings are set and after addTitleHead is called
    param = 'view'
    values = ["0.20", str(self.bottom), "0.90", str(self.top)]
    value = ', '.join(values)
    self.addBulkParams([param],[value])
   

  def addAxes(self):
    axes = self.axes
    
    params = ['xaxis label', 'yaxis label']
    values = [axes['label'][0], axes['label'][1]]
    
    self.addBulkParams(params,values)

  def stringsToParams(self):
    # TODO WILL NEED TO BE EDITTED !!!!
    for string in self.strings:
      params = ['WITH STRING', 'STRING LOCTYPE', 'STRING CHAR SIZE', 'STRING', 'STRING DEF', 'STRING']
      values = ["", 'VIEW', '0.9', '0.2, ' + str(string['position']), string['string'], 'on']
      self.addBulkParams(params,values)

  def addData(self):
    i = 0;
    scatter = 0;
    for dataSet in self.dataSets:
      paramPrefix = "s" +str(i)
      if 'color' in dataSet:
        color = str(dataSet['color'])
      else:
        color = str(i+1)
      params = ['type']
      if dataSet['plot'] == 'scatter':
        params += ['line type', 'symbol','symbol color', 'symbol pattern', 'symbol fill color','symbol fill pattern','errorbar color','symbol size']
        values = ['xydy','0', str(scatter+1), color, '2', color,'1',color,'.75']
        scatter += 1;
      if dataSet['plot'] == 'line':
        params += ['line type','line color','symbol']
        values = ['xy','1',color,'0']
      i += 1;
      params = [(paramPrefix + " " + param) for param in params]
      
      self.addBulkParams(params, values)


  def write(self):
    self.addTitleHead()
    self.addViews()
    self.addAxes()
    self.stringsToParams()
    self.addData()
#    self.writeOut('test.xmgr')

  def writeOut(self,file):
    target = open(file,'w');
    for line in self.parameters:
      target.write(line)
    i = 0;
    for dataSet in self.dataSets:
      if i != 0:
        target.write('&\n')
      string = "@type ";
      if dataSet['plot'] == 'scatter':
        string += "xydy"
      else:
         string += "xy"
      string +="\n"
      target.write(string) 
      for data in dataSet['data']:
        data += "\n"
        target.write(data)
        i+=1
    target.close()
  

# This is not needed once data is brought in from java
def createData(dataPts,type="linear", m=-3, b=6):
  #type = linear, log
  data = []
  
  if type=="linear":
    for i in range(dataPts):
      x = i/10.0;
      y = m*x+b+random.gauss(0,.5)
      stderr = random.gauss(0,.3)
      datum = [x,y,stderr]
      datum = [str(i) for i in datum]
      line = " ".join(datum)
      data.append(line)
  return data

def test():
    print "Hello World"  

#writer = Writer()

#fakeData1 = {'data' : createData(100,m=2,b=0), 'plot': 'scatter'}
#fakeData2 = {'data' : createData(100,m=1.0, b= .5), 'plot': 'scatter'}
#writer.addDataSets([fakeData1,fakeData2])
#writer.write()
#writer.writeOut('./test.xmgr')
      
#def testRead(list):
#    for i in list:
#        print i


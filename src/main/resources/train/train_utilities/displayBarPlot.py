# The purpose of this script is to use the information passed by analyze_performance.py to generate some bar plots that'll
# demonstrate the accuracy from the current guesser and the neural network guesser.


import sys
import os
import matplotlib.pyplot as plt                                                                                                                                         
import numpy as np

def displayBarPlot():
    mode = sys.argv[1]
    nFields = int(sys.argv[2])
    rmseList = [[float(val) for val in pair.split(":") if val] for pair in sys.argv[3].split(',') if pair]

    mode = mode.upper()
    networkGuess = []
    iterativeGuess = []
    nTests = 15 if len(rmseList) > 20 else len(rmseList)
    for i in range(nTests):
        iterativeRMSE, networkRMSE = rmseList[i]
        iterativeGuess.append(iterativeRMSE)
        networkGuess.append(networkRMSE)

    ind = np.arange(nTests)
    width = 0.35
    plt.bar(ind,iterativeGuess, width, label='Old Guesser', color='r')
    plt.bar(ind + width,networkGuess, width, label='Neural Network Guesser', color='k')
    
    plt.xlabel('Tests')
    plt.ylabel('RMSE')
    pluralField = 'field' if nFields == 1 else 'fields'
    plt.title('{} RMSE Scores'.format(mode))

    plt.xticks(ind + width / 2, [str(i) for i in range(1,nTests+1)])
    plt.legend(loc='best')
    plt.savefig("./{0}_{1}Fields.png".format(mode,str(nFields)))
    plt.show()

displayBarPlot()

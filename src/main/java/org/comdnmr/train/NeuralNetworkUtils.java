package org.comdnmr.fit.train;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;
import org.ojalgo.matrix.store.MatrixStore;
import org.ojalgo.ann.ArtificialNeuralNetwork;
import org.ojalgo.ann.NetworkBuilder;
import org.ojalgo.structure.Access1D;
import java.util.HashMap;

/**
 * Purpose - This class is used to contain utilities that could be of useful for
 * training, testing, saving an artificial neural network.
 *
 *
 * @author teddycolon
 */
public class NeuralNetworkUtils {

    public static String saveNeuralNetwork(String saveDir, ArtificialNeuralNetwork trainedNetwork, HashMap networkInfo) throws IOException {
        String name = (String) networkInfo.get("name");
        String SaveNeuralNetworkFile = saveDir + "/Saved_ANN_{}.txt".replace("{}", name);
        FileWriter fileWriter = new FileWriter(SaveNeuralNetworkFile);
        PrintWriter printWriter = new PrintWriter(fileWriter);
        int[] nodesPerLayer = (int[]) networkInfo.get("nodesPerLayer");
        ArtificialNeuralNetwork.Activator[] activators = (ArtificialNeuralNetwork.Activator[]) networkInfo.get("activators");
        printWriter.printf("Neurons Per Layer: %s", Arrays.toString(nodesPerLayer));
        for (int layer = 0, limit = (int) networkInfo.get("nLayers"); layer < limit; layer++) {
            int inNodes = nodesPerLayer[layer];
            int outNodes = nodesPerLayer[layer + 1];
            // weights
            printWriter.format("\nLayer (%d)\n", layer);
            printWriter.println("Weights:");
            for (int input = 0; input < inNodes; input++) {
                for (int output = 0; output < outNodes; output++) {
                    double weight = trainedNetwork.getWeight(layer, input, output);
                    printWriter.printf("%-3d %-3d %3d %21f\n", layer, input, output, weight);
                }
            }
            // biases
            printWriter.println("");
            printWriter.println("Biases:");
            for (int output = 0; output < outNodes; output++) {
                double bias = trainedNetwork.getBias(layer, output);
                printWriter.printf("%-3d %-3d %15f\n", layer, output, bias);
            }

        }
        printWriter.println("");
        printWriter.printf("Activators: %s\n", Arrays.toString(activators));
        if (networkInfo.containsKey("scaleValues")) {
            printWriter.printf("Scale values = %s", networkInfo.get("scaleValues"));
        }

        printWriter.close();
        return SaveNeuralNetworkFile;
    }

    public static int[] getExampleLineInfo(String exampleLine) {
        String[] colonSplit = exampleLine.split(":");
        int[] exampleLineInfo = new int[colonSplit.length];
        int ind = 0;
        for (String str : colonSplit) {
            String[] commaSplit = str.split(",");
            exampleLineInfo[ind] = commaSplit.length;
            ind++;
        }
        return exampleLineInfo;
    }

    public static List convertExampleLine(String fileLine, int nInputParams, int nOutputParams) throws Exception {
        List exampleLine = new ArrayList(2);
        List<Double> inputValues = new ArrayList(nInputParams);
        List<Double> targetValues = new ArrayList(nOutputParams);
        double convertedValue;

        String[] arrStr = fileLine.split(":");
        for (String arr : arrStr) {
            String[] arrOfStringNums = arr.split(",");
            for (int i = 0, nNums = arrOfStringNums.length; i < nNums; i++) {
                convertedValue = Double.parseDouble(arrOfStringNums[i]);
                if (nNums == nOutputParams) {
                    targetValues.add(convertedValue);
                } else if (nNums == nInputParams) {
                    inputValues.add(convertedValue);
                } else {
                    throw new Exception("Number of values for example line ({}) is invalid.".replace("{}", fileLine));
                }
            }
        }
        exampleLine.add(inputValues);
        exampleLine.add(targetValues);
        return exampleLine;
    }

    public static double errorCalculation(MatrixStore<Double> networkOutput, List<Double> expectedOutput) {
        double error;
        double sumSq = 0.0;
        for (int i = 0; i < networkOutput.size(); i++) {
            error = expectedOutput.get(i) - networkOutput.get(i);
            sumSq += (error * error);
        }
        return sumSq;
    }

    public static MyResult trainNetwork(List trainFile, List validFile, NetworkBuilder neuralBuilder) throws Exception {
        
        String line;
        List trainLineAsList; // List<Double> or List<List<Double>>
        List<Double> trainInput,
                targetOutput;
        double bestTrainingError = 0.0001;
        double currMSE = 0.0;
        int epoch = 0,
                maxEpochs = 300;
        Access1D<Double> input1D,
                target1D;

        do { // iterate over training set
            
            // training
            int i=0;
            while (i < trainFile.size()) { // training set 
                line = (String) trainFile.get(i);

                // grab test line in a processable format
                int[] exampleLineInfo = getExampleLineInfo(line);
                int nInputNodes = exampleLineInfo[0];
                int nOutputNodes = exampleLineInfo[1];
                trainLineAsList = convertExampleLine(line, nInputNodes, nOutputNodes);

                trainInput = (List<Double>) trainLineAsList.get(0);
                targetOutput = (List<Double>) trainLineAsList.get(1);

                // conversion for simplicity to Access1D
                input1D = Access1D.wrap(trainInput);
                target1D = Access1D.wrap(targetOutput);
                neuralBuilder.train(input1D, target1D); // supervised training function
                i++;
            }

            // validation
            if (validFile != null) {
                int j=0;
                int nGuessedValues = 0;
                int nExamples = 0;
                while (j < validFile.size()) {
                    ArtificialNeuralNetwork trainedNetwork = neuralBuilder.get();
                    line = (String) validFile.get(j);
                    int[] exampleLineInfo = getExampleLineInfo(line);
                    int nInputNodes = exampleLineInfo[0];
                    int nOutputNodes = exampleLineInfo[1];
                    List validationLineAsList = convertExampleLine(line, nInputNodes, nOutputNodes);
                    List<Double> validationInput = (List<Double>) validationLineAsList.get(0);
                    List<Double> validationTarget = (List<Double>) validationLineAsList.get(1);
                    Access1D validWrapInput = Access1D.wrap(validationInput);
                    MatrixStore validNetworkOuput = trainedNetwork.invoke(validWrapInput);
                    nGuessedValues += validationTarget.size();
                    nExamples += 1;
                    currMSE += errorCalculation(validNetworkOuput, validationTarget) / nGuessedValues;
                    j++;
                }   
                currMSE /= nExamples;
            }
            epoch++;
            //System.out.println("mse : " + currMSE);
        } while (epoch < maxEpochs && currMSE > bestTrainingError);
        
//        System.out.println("1.5_ current MSE : " + currMSE);
        String networkString = neuralBuilder.toString();

        return new MyResult<>(neuralBuilder, networkString, currMSE);
    }

    public static double testNetwork(Scanner testFile, ArtificialNeuralNetwork trainedNetwork, int nInputNeurons, int nOutputNeurons) throws Exception {
        String line;
        List testLineAsList;
        List<Double> testInput;
        List<Double> expectedTestOutput;
        Access1D testInput1D;
        MatrixStore actualOutput;
        double sumSqTest = 0.0;
        int nTest = 0;
        double testingMeanSqError = 0.0;

        // Input and target output are expected to be normalized
        while (testFile.hasNextLine()) {
            // GET LINE AND GRAB INPUT, EXPECTED OUTPUT, SCALE VALUE
            line = testFile.nextLine();
            int[] exampleLineInfo = getExampleLineInfo(line);
            int nInputNodes = exampleLineInfo[0];
            int nOutputNodes = exampleLineInfo[1];
            testLineAsList = convertExampleLine(line, nInputNodes, nOutputNodes);
            testInput = (List<Double>) testLineAsList.get(0);
            expectedTestOutput = (List<Double>) testLineAsList.get(1);

            // WRAP INPUT IN Access1D FOR NETWORK INPUT
            testInput1D = Access1D.wrap(testInput);

            // INVOKE NETWORK AND COMPUTE ERROR
            actualOutput = trainedNetwork.invoke(testInput1D);
//            BasicLogger.debug("{} <> {}", revertScaleOutput(actualOutput, paramScaleVals), expectedTestOutput.toString());

            sumSqTest += errorCalculation(actualOutput, expectedTestOutput);
            nTest += expectedTestOutput.size();
        }
        testingMeanSqError = sumSqTest / nTest;
//        System.out.println("MSE(test) : " + testingMeanSqError);
//        System.out.println("RMSE(test) : " + Math.sqrt(testingMeanSqError));
        return testingMeanSqError; // MSE(test)
    }

    public static double[] getSavedNetworkOutput(String savedNetworkFile, List yValues) throws Exception {
        ANNLoader ANN = new ANNLoader(savedNetworkFile); // cpmg fast, slow for 1 and 2 fields, cestr1rhoperturbation
        System.out.println("-- Loaded saved neural network!");
        ArtificialNeuralNetwork trainedNetwork = ANN.getTrainedNetwork();
        Access1D wrappedInput = Access1D.wrap(yValues); // input vals already scaled
        MatrixStore<Double> result = trainedNetwork.invoke(wrappedInput); //scaled output pars; need to revert
        System.out.println("-- Network invoked!");
        return result.toRawCopy1D();
    }

}

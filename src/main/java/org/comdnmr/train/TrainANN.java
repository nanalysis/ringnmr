/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.fit.train;

import static org.comdnmr.fit.train.NeuralNetworkUtils.trainNetwork;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Scanner;
import java.util.Random;
import org.ojalgo.ann.ArtificialNeuralNetwork;
import static org.ojalgo.ann.ArtificialNeuralNetwork.Activator.IDENTITY;
import static org.ojalgo.ann.ArtificialNeuralNetwork.Activator.RECTIFIER;
import static org.ojalgo.ann.ArtificialNeuralNetwork.Error.HALF_SQUARED_DIFFERENCE;
import org.ojalgo.ann.NetworkBuilder;

/**
 *
 * @author teddycolon
 */
public abstract class TrainANN {

    private String tFilePathString = "";
    private String vFilePathString = "";
    private String scaleValues = "";
    private final List<String> arrTrainFile = new ArrayList();
    private final List<String> arrValidFile = new ArrayList();

    public TrainANN(String trainFileName) throws FileNotFoundException {
        // tFilePathString should be an absolute path string to file
        this.tFilePathString = trainFileName;
        File trainFile = new File(tFilePathString);
        if (trainFile.exists()) {
            Scanner scanTrainFile = new Scanner(trainFile);
            while (scanTrainFile.hasNextLine()) {
                this.arrTrainFile.add(scanTrainFile.nextLine());
            }
        }
    }

    // constructor that accepts validation file name
    public TrainANN(String trainFileName, String validationFileName) throws FileNotFoundException {
        this.tFilePathString = trainFileName;
        this.vFilePathString = validationFileName;
        File trainFile = new File(tFilePathString);
        File validFile = new File(vFilePathString);
        if (trainFile.exists() && validFile.exists()) {
            Scanner scanTrainFile = new Scanner(trainFile);
            Scanner scanValidFile = new Scanner(validFile);
            while (scanTrainFile.hasNextLine()) this.arrTrainFile.add(scanTrainFile.nextLine());
            while (scanValidFile.hasNextLine()) this.arrValidFile.add(scanValidFile.nextLine());
        }
    }

    /**
     * Purpose - Setter method to allow user to write scale values into the
     * saved network file.
     *
     * @param values colon delimited string which contains parameter name, max
     * and min values for each parameter delimited by semicolons and commas,
     * respectively.
     *
     * For example:
     *
     * Parameters = ["kex", "Rex"], Scale values = [(max1, min1), (max2,min2)]
     * then scale values string is "Kex;max1,min1:Rex;max2,min2"
     */
    public void setScaleValues(String values) {
        this.scaleValues = values;
    }

    private MyResult trainExecutor(String trainLabel, int[] nodesPerLayer, ArtificialNeuralNetwork.Activator[] activators, double learningRate) throws FileNotFoundException, Exception {

        int nInputNeurons = nodesPerLayer[0];
        NetworkBuilder builder = null, trainedBuilder = null;

        // INITIALIZING NETWORK DESIGN
        builder = ArtificialNeuralNetwork.builder(nInputNeurons, Arrays.copyOfRange(nodesPerLayer, 1, nodesPerLayer.length));
        builder.activators(activators).error(HALF_SQUARED_DIFFERENCE).rate(learningRate);

        // TRAINING NETWORK
        MyResult trainingResult = trainNetwork(arrTrainFile, arrValidFile, builder);
        trainedBuilder = (NetworkBuilder) trainingResult.getFirst();

        // GRABBING TRAINED NETWORK
        ArtificialNeuralNetwork trainedNetwork = trainedBuilder.get();
        
        // STORING NETWORK INFO
        HashMap networkInfo = new HashMap<>(5);
        networkInfo.put("name", trainLabel);
        networkInfo.put("nodesPerLayer", nodesPerLayer);
        networkInfo.put("nLayers", trainedNetwork.countCalculationLayers());
        networkInfo.put("activators", activators);
        if (!scaleValues.isEmpty()) {
            networkInfo.put("scaleValues", scaleValues);
        }

        double trainingMSE = (double) trainingResult.getThird();

        return new MyResult<>(trainedNetwork, networkInfo, trainingMSE);
    }

    public MyResult trainingRun(String trainLabel, int nInputNeurons, int nOutputNeurons) throws Exception {
        int nCollectiveLayers = 5;
        ArtificialNeuralNetwork.Activator[] activators = {RECTIFIER, RECTIFIER, RECTIFIER, IDENTITY};
        int[] nodesForAllLayers = new int[nCollectiveLayers];
        nodesForAllLayers[0] = nInputNeurons;
        nodesForAllLayers[nCollectiveLayers - 1] = nOutputNeurons;
        Random randGen = new Random();
        int maxTrainModels = 5,
                currTrainModel = 0;
        double initLearningRate = 0.2;
        MyResult bestResult = null;
        System.out.println("0_ train file : " + tFilePathString);
        System.out.println("0_ validation file : " + vFilePathString);
        System.out.println("0_ nInputNeurons >>> " + nInputNeurons);
        System.out.println("0_ nOutputNeurons >>> " + nOutputNeurons);
        System.out.println("------------------------------------------------------------");

        while (currTrainModel < maxTrainModels) {

            // initializing the number of nodes in the hidden layers
            for (int hL = 1, lastLayer = nCollectiveLayers - 1; hL < lastLayer; hL++) {
                int diff = nodesForAllLayers[hL-1] - nOutputNeurons;
                int upperBound = (diff <= 0) ?  nOutputNeurons : diff;
                int randInt = randGen.nextInt(upperBound) + nOutputNeurons;
                
                // EDIT: Change so that the number of nodes per hidden layer decrease progressively
                nodesForAllLayers[hL] = randInt;
            }
            
            System.out.println("1_ nodes for all layers >>> " + Arrays.toString(nodesForAllLayers));
            int[] nodes = (int[]) nodesForAllLayers.clone();
            MyResult result = trainExecutor(trainLabel, nodes, activators, initLearningRate);
            
            // running training model with different learning rates to find which learning rate produces smallest error
            for (double rate = initLearningRate + 0.01, limit = 0.3; rate < limit; rate += 0.01) {
                MyResult tempResult = trainExecutor(trainLabel, nodes, activators, rate);
                if ((double) tempResult.getThird() < (double) result.getThird()) {
                    result = new MyResult<>((ArtificialNeuralNetwork) tempResult.getFirst(), (HashMap) tempResult.getSecond(), (double) tempResult.getThird());
                }
            }

            // update best result based on the result with the lowest error
            double resultScore = (double) result.getThird();
            double bestResultScore = (bestResult == null) ? resultScore : (double) bestResult.getThird();
            
            if ((bestResult == null) || (resultScore < bestResultScore)) {
//                bestResult = result.cloning();
                bestResult = new MyResult<>((ArtificialNeuralNetwork) result.getFirst(), (HashMap) result.getSecond(), (double) result.getThird());
//                System.out.println("------------------------------------------------------------");
            }
            System.out.println("__3.1_ result score : " + resultScore);
            System.out.println("__3.2_ best result score : " + bestResultScore);
//            String printStr = "{1} nodes per layer: ".replace("{1}", String.valueOf(currTrainModel));
//            System.out.println(printStr + Arrays.toString(nodesForAllLayers));
            currTrainModel++;
        }
        System.out.println("------------------------------------------------------------");
//        ArtificialNeuralNetwork trainedNetwork = (ArtificialNeuralNetwork) bestResult.getFirst();
        HashMap networkInfo = (HashMap) bestResult.getSecond();
        int[] nPL = (int[]) networkInfo.get("nodesPerLayer");
        System.out.println("4_ Chosen BestResult nodes per layer >>> " + Arrays.toString(nPL));
        System.out.println("------------------------------------------------------------");

        return bestResult;
    }

}

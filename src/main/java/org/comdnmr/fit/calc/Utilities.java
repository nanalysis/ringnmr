
package org.comdnmr.fit.calc;

/**
 *
 * @author brucejohnson
 */
public class Utilities {
    
    public static double[][] copy2DArray(double[][] inData) {
        double[][] outData = new double[inData.length][];
        for (int i=0;i<inData.length;i++) {
            outData[i] = inData[i].clone();
        }
        return outData;        
    }
}

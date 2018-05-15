package org.comdnmr.cpmgfit2.calc;

import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author Bruce Johnson
 */
public class ResidueData {

    ExperimentData expData;
    String resNum;
    double[][] xValues;
    double[] errValues;
    double[] yValues;
    String[] peakRefs;

    public ResidueData(ExperimentData expData, String residueNum, double[][] x, double[] y, double[] err) {
        this.resNum = residueNum;
        this.expData = expData;
        this.xValues = DataUtil.clone2DArray(x);
        this.yValues = y.clone();
        this.errValues = err.clone();
    }

    public ResidueData(ExperimentData expData, String residueNum, double[][] x, double[] y, double[] err, String[] peakRefs) {
        this.resNum = residueNum;
        this.expData = expData;
        this.xValues = DataUtil.clone2DArray(x);
        this.yValues = y.clone();
        this.errValues = err.clone();
        this.peakRefs = peakRefs.clone();
    }

    public ResidueData(ExperimentData expData, String residueNum, List<Double> xValueList, List<Double> yValueList, List<Double> errValueList) {
        this.expData = expData;
        this.resNum = residueNum;
        int nValues = xValueList.size();
        this.xValues = new double[1][nValues];
        this.yValues = new double[nValues];
        this.errValues = new double[nValues];
        for (int i = 0; i < nValues; i++) {
            xValues[0][i] = xValueList.get(i);
            yValues[i] = yValueList.get(i);
            errValues[i] = errValueList.get(i);
        }
    }

    public ResidueData(ExperimentData expData, String residueNum, List<Double>[] xValueList, List<Double> yValueList, List<Double> errValueList) {
        this.expData = expData;
        this.resNum = residueNum;
        int nValues = yValueList.size();
        int nX = xValueList.length;
        this.xValues = new double[nX][nValues];
        this.yValues = new double[nValues];
        this.errValues = new double[nValues];
        for (int i = 0; i < nValues; i++) {
            for (int j = 0; j < nX; j++) {
                xValues[j][i] = xValueList[j].get(i);
            }
            yValues[i] = yValueList.get(i);
            errValues[i] = errValueList.get(i);
        }
    }

    public ResidueData(ExperimentData expData, String residueNum, List<Double> xValueList, List<Double> yValueList, List<Double> errValueList, List<String> peakRefList) {
        this.expData = expData;
        this.resNum = residueNum;
        int nValues = xValueList.size();
        this.xValues = new double[1][nValues];
        this.yValues = new double[nValues];
        this.errValues = new double[nValues];
        this.peakRefs = new String[nValues];
        for (int i = 0; i < nValues; i++) {
            xValues[0][i] = xValueList.get(i);
            yValues[i] = yValueList.get(i);
            errValues[i] = errValueList.get(i);
            peakRefs[i] = peakRefList.get(i);
        }
    }

    public double[][] getXValues() {
        return xValues;
    }

    public double[] getYValues() {
        return yValues;
    }

    public double[] getErrValues() {
        return errValues;
    }

    public class DataValue {

        int index;
        ResidueData resInfo;

        public DataValue(ResidueData resInfo, int index) {
            this.resInfo = resInfo;
            this.index = index;
        }

        public int getIndex() {
            return index;
        }

        public double getX0() {
            return resInfo.xValues[0][index];
        }
        
        public double getX1() {
            double x1 = resInfo.xValues[0][index];
            if (resInfo.xValues.length > 1) {
                x1 = resInfo.xValues[1][index];
            } 
            return x1;
        }

        public double getY() {
            return resInfo.yValues[index];
        }

        public double getError() {
            return resInfo.errValues[index];
        }

        public String getPeak() {
            String peak = "";
            if (resInfo.peakRefs != null) {
                peak = resInfo.peakRefs[index];
            }
            return peak;
        }

        public String getName() {
            return resInfo.getDataName();
        }

        public String getResidue() {
            return resInfo.resNum;
        }
    }

    public ArrayList<DataValue> getDataValues() {
        ArrayList<DataValue> dataValues = new ArrayList<>();
        for (int i = 0; i < xValues[0].length; i++) {
            dataValues.add(new DataValue(this, i));
        }
        return dataValues;
    }

    public String getDataName() {
        return expData.getName();
    }

    public ExperimentData getExperimentData() {
        return expData;
    }

}

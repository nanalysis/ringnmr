package org.comdnmr.cpmgfit2.calc;

import java.util.ArrayList;

/**
 *
 * @author Bruce Johnson
 */
public class ResidueData {

    int resNum;
    double[] xValues;
    double[] errValues;
    double[] yValues;
    String[] peakRefs;

    public ResidueData(double[] x, double[] y, double[] err) {
        this.xValues = x.clone();
        this.yValues = y.clone();
        this.errValues = err.clone();
    }

    public double[] getXValues() {
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

        public double getX() {
            return resInfo.xValues[index];
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
    }

    public ArrayList<DataValue> getDataValues() {
        ArrayList<DataValue> dataValues = new ArrayList<>();
        for (int i = 0; i < xValues.length; i++) {
            dataValues.add(new DataValue(this, i));
        }
        return dataValues;
    }

}

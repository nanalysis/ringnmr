/*
 * CoMD/NMR Software : A Program for Analyzing NMR Dynamics Data
 * Copyright (C) 2018-2019 Bruce A Johnson
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package org.comdnmr.data;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;
import org.comdnmr.util.DataUtil;
import org.nmrfx.chemistry.Atom;
import org.nmrfx.chemistry.Entity;
import org.nmrfx.chemistry.MoleculeBase;
import org.nmrfx.chemistry.MoleculeFactory;
import org.nmrfx.chemistry.RelaxationData;
import org.nmrfx.chemistry.RelaxationData.relaxTypes;
import org.nmrfx.chemistry.Residue;

/**
 *
 * @author Bruce Johnson
 */
public class ExperimentalData {

    Experiment experiment;
    DynamicsSource dynSource;
    RelaxationData relaxData;
    double[][] xValues;
    double[] errValues;
    double[] yValues;

    public ExperimentalData(Experiment experiment, DynamicsSource dynSource,
            double[][] x, double[] y, double[] err) {
        this.experiment = experiment;
        this.dynSource = dynSource;
        this.xValues = DataUtil.clone2DArray(x);
        this.yValues = y.clone();
        this.errValues = err.clone();
    }

    public ExperimentalData(Experiment experiment, DynamicsSource dynSource,
            List<Double> xValueList, List<Double> yValueList, List<Double> errValueList) {
        this.experiment = experiment;
        this.dynSource = dynSource;
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

    public ExperimentalData(Experiment experiment, DynamicsSource dynSource,
            List<Double>[] xValueList, List<Double> yValueList,
            List<Double> errValueList) {
        this.experiment = experiment;
        this.dynSource = dynSource;
        int nValues = yValueList.size();
        int nX = xValueList.length;
//        System.out.println("make res data " + nX);
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

    public ExperimentalData(Experiment expData, RelaxationData relaxData) {
        this.experiment = expData;
        this.relaxData = relaxData;
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
        ExperimentalData resInfo;

        public DataValue(ExperimentalData resInfo, int index) {
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

        public double getX2() {
            double x2 = resInfo.xValues[0][index];
            if (resInfo.xValues.length > 2) {
                x2 = resInfo.xValues[2][index];
            }
            return x2;
        }

        public double getY() {
            return resInfo.yValues[index];
        }

        public double getError() {
            return resInfo.errValues[index];
        }

        public String getRefPeak() {
            return resInfo.dynSource.peak.getName();
        }

        public int getPeak() {
            return resInfo.dynSource.peak.getIdNum();
        }

        public String getName() {
            return resInfo.getDataName();
        }

        public String getResidue() {
            Entity entity = resInfo.dynSource.atoms[0].getEntity();
            if (entity instanceof Residue) {
                Residue residue = (Residue) entity;
                return residue.getNumber();
            } else {
                return "1";
            }
        }

        public String getResName() {
            Entity entity = resInfo.dynSource.atoms[0].getEntity();
            if (entity instanceof Residue) {
                Residue residue = (Residue) entity;
                return residue.getName();
            } else if (entity != null) {
                return entity.getName();
            } else {
                return "";
            }
        }
    }

    public ArrayList<DataValue> getDataValues() {
        ArrayList<DataValue> dataValues = new ArrayList<>();
        if (xValues != null) {
            for (int i = 0; i < xValues[0].length; i++) {
                dataValues.add(new DataValue(this, i));
            }
        }
        return dataValues;
    }

    public String getDataName() {
        return experiment.getName();
    }

    public Experiment getExperimentData() {
        return experiment;
    }

    public DynamicsSource getDynSource() {
        return dynSource;
    }

}

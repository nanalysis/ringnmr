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
import java.util.Map;
import java.util.Set;
import org.comdnmr.util.DataUtil;
import org.nmrfx.chemistry.Atom;
import org.nmrfx.chemistry.Entity;
import org.nmrfx.chemistry.MoleculeBase;
import org.nmrfx.chemistry.MoleculeFactory;
import org.nmrfx.chemistry.Residue;

/**
 *
 * @author Bruce Johnson
 */
public class ResidueData {

    ExperimentData expData;
    DynamicsSource dynSource;
    double[][] xValues;
    double[] errValues;
    double[] yValues;

    public ResidueData(ExperimentData expData, DynamicsSource dynSource,
            double[][] x, double[] y, double[] err) {
        this.expData = expData;
        this.dynSource = dynSource;
        this.xValues = DataUtil.clone2DArray(x);
        this.yValues = y.clone();
        this.errValues = err.clone();
    }

    public ResidueData(ExperimentData expData, DynamicsSource dynSource,
            List<Double> xValueList, List<Double> yValueList, List<Double> errValueList) {
        this.expData = expData;
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

    public ResidueData(ExperimentData expData, DynamicsSource dynSource,
            List<Double>[] xValueList, List<Double> yValueList,
            List<Double> errValueList) {
        this.expData = expData;
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
            } else {
                return entity.getName();
            }
        }
        
        public String getT1() {
            MoleculeBase mol = MoleculeFactory.getActive();
            mol.getAtomArray();
            Atom[] atoms = resInfo.dynSource.atoms; 
            String t1 = "-";
            for (Atom atom : atoms) {
                String aName = atom.getFullName();
                if (aName.contains(resInfo.expData.nucleusName) && mol.getAtom(aName) != null) {
                    Set<String> fields = (Set<String>) mol.getAtom(aName).getProperty("T1fields");
                    if (fields != null) {
                        Map<String, String> expVals = (Map<String, String>) mol.getAtom(aName).getProperty("T1" + fields.toArray()[0] + "_1" + "results"); 
                        if (expVals != null) {
                            t1 = expVals.get("T1Val");
                        }
                    }
                }
            }
            return t1;
        }

        public String getT2() {
            MoleculeBase mol = MoleculeFactory.getActive();
            mol.getAtomArray();
            Atom[] atoms = resInfo.dynSource.atoms; 
            String t2 = "-";
            for (Atom atom : atoms) {
                String aName = atom.getFullName();
                if (aName.contains(resInfo.expData.nucleusName) && mol.getAtom(aName) != null) {
                    Set<String> fields = (Set<String>) mol.getAtom(aName).getProperty("T2fields");
                    if (fields != null) {
                        Map<String, String> expVals = (Map<String, String>) mol.getAtom(aName).getProperty("T2" + fields.toArray()[0] + "_1" + "results"); 
                        if (expVals != null) {
                            t2 = expVals.get("T2Val");
                        }
                    }
                }
            }
            return t2;
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
    
    public DynamicsSource getDynSource() {
        return dynSource;
    }

}

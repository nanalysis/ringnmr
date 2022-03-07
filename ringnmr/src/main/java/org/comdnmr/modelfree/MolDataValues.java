/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.modelfree;

import org.comdnmr.data.DynamicsSource;
import org.comdnmr.modelfree.models.MFModel;
import org.nmrfx.chemistry.Atom;

import java.util.ArrayList;
import java.util.List;
import java.util.Optional;

/**
 *
 * @author brucejohnson
 */
public class MolDataValues {

    Atom atom;
    final String specifier;
    final double[] vector = new double[3];
    final List<RelaxDataValue> dataValues = new ArrayList<>();
    MFModel model;
    double[][] jValues = null;

    public MolDataValues(String specifier, double[] vector, DynamicsSource dynSourceFactory) {
        this.specifier = specifier;
        Optional<Atom> atomOpt = dynSourceFactory.getAtom(null, specifier, 0, false);
        atom = atomOpt.orElse(null);
        System.arraycopy(vector, 0, this.vector, 0, 3);
    }

    public MolDataValues(Atom atom, double[] vector) {
        this.atom = atom;
        this.specifier = atom.getFullName();
        System.arraycopy(vector, 0, this.vector, 0, 3);
    }

    public MolDataValues(Atom atom) {
        this.atom = atom;
        this.specifier = atom.getFullName();
    }

    @Override
    public String toString() {
        StringBuilder sBuilder = new StringBuilder();
        for (var dataValue : dataValues) {
            sBuilder.append(dataValue.toString());
            sBuilder.append("\n");
        }
        return sBuilder.toString();
    }

    public void addData(RelaxDataValue value) {
        dataValues.add(value);
    }

    public List<RelaxDataValue> getData() {
        return dataValues;
    }

    public void setTestModel(MFModel model) {
        this.model = model;
    }

    public MFModel getTestModel() {
        return model;
    }

    public double[] getVector() {
        return vector;
    }

    public double[][] calcJ() {
        int nDataValues = dataValues.size();
        int nFreq = 1 + 2 * nDataValues;
        double[][] result = new double[3][nFreq];
        int iField = 0;
        for (RelaxDataValue relaxDataValue: dataValues) {
            RelaxEquations relaxEq = relaxDataValue.relaxObj;
            double field = relaxEq.getSF();
            double r1 = relaxDataValue.R1;
            double r2 = relaxDataValue.R2;
            double noe = relaxDataValue.NOE;
            double r1Err = relaxDataValue.R1err;
            double r2Err = relaxDataValue.R2err;
            double noeErr = relaxDataValue.NOEerr;
            double sigma = (noe - 1.0) * r1 * RelaxEquations.GAMMA_N / RelaxEquations.GAMMA_H;
            double sigmaErr = sigma * Math.sqrt(Math.pow((noeErr / (noe - 1.0)), 2) + Math.pow((r1Err / r1), 2));

            double d2 = relaxEq.getD2();
            double c2 = relaxEq.getC2();

            double j87H = 4.0 * sigma / (5.0 * d2);
            double j87Herr = 4.0 * sigmaErr / (5.0 * d2);

            double jNMul = 4.0 / (3.0 * d2 + 4.0 * c2);
            double jN = (r1 - 1.249 * sigma) * jNMul;
            double jNerr = jNMul * Math.sqrt(Math.pow(r1Err, 2) + Math.pow(1.249 * sigmaErr, 2));

            double j0Mul = 6.0 / (3.0 * d2 + 4.0 * c2);
            double j0 = j0Mul * (r2 - 0.5 * r1 - 0.454 * sigma);
            double j0Err = j0Mul * Math.sqrt(Math.pow(r2Err, 2) + Math.pow(0.5 * r1Err, 2) + Math.pow(0.454 * sigmaErr, 2));
            result[0][0] += j0;
            result[1][0] = 0.0;
            result[2][0] += j0Err *j0Err;

            result[0][iField*2+1] = j87H;
            result[1][iField*2+1] = 0.87 * relaxEq.getWI();
            result[2][iField*2+1] = j87Herr;

            result[0][iField*2+2] = jN;
            result[1][iField*2+2] = 0.87 * relaxEq.getWS();
            result[2][iField*2+2] = jNerr;
            iField++;
        }
        result[0][0] /= nDataValues;
        result[2][0] = Math.sqrt(result[2][0]);
        return result;
    }

    double[][] getJValues() {
        if (jValues == null) {
            jValues = calcJ();
        }
        return jValues;
    }
}

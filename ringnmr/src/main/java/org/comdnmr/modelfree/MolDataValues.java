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
        var dataOpt = dataValues.stream().findFirst();
        if (dataOpt.isPresent()) {
            if (dataOpt.get() instanceof DeuteriumDataValue) {
                return SpectralDensityCalculator.calcJDeuterium(dataValues);
            } else {
                return SpectralDensityCalculator.calcJR1R2NOE(dataValues);
            }
        }
        return new double[0][0];
    }

    double[][] getJValues() {
        if (jValues == null) {
            jValues = calcJ();
        }
        return jValues;
    }
}

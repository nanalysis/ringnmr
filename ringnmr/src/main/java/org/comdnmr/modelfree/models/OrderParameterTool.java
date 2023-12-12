package org.comdnmr.modelfree.models;

import org.comdnmr.modelfree.FitDeuteriumModel;
import org.comdnmr.modelfree.MolDataValues;
import org.nmrfx.chemistry.Atom;
import org.nmrfx.chemistry.MoleculeBase;
import org.nmrfx.chemistry.MoleculeFactory;
import org.nmrfx.chemistry.relax.*;

import java.util.Collections;
import java.util.Map;

public class OrderParameterTool {

    private OrderParameterTool() {

    }

    public static void calculateOrderParameters() {
        var data = FitDeuteriumModel.getData(false);
        data.entrySet().stream().sorted(Map.Entry.comparingByKey()).forEach(e -> {
            MolDataValues resData = e.getValue();
            Atom atom = resData.getAtom();
            String key = e.getKey();
            resData.clearJ();
            double[][] jValues = resData.getJValues();
            SpectralDensity spectralDensity = new SpectralDensity(key, jValues);
            atom.addSpectralDensity(key, spectralDensity);
        });
    }
    public static void interpolateRates(double field2) {
        MoleculeBase moleculeBase = MoleculeFactory.getActive();
        var relaxMap = moleculeBase.relaxationSetMap();
        Double lower = null;
        Double upper = null;
        for (var entry : relaxMap.entrySet()) {
            double field = entry.getValue().field();
            if ((field2 - field) > 1.0) {
                if ((lower == null) || (field > lower)) {
                    lower = field;
                }
            } else if ((field2 - field) < -1.0) {
                if ((upper == null) || (field < upper)) {
                    upper = field;
                }
            }
        }
        System.out.println("interpolate " + lower + " " + field2 + " " + upper);
        if ((lower != null) && (upper != null)) {
            interpolateRates(lower, field2, upper);
        }
    }

    public static void interpolateRates(double field1, double field2, double field3) {
        var data = FitDeuteriumModel.getData(false);
        data.entrySet().stream().sorted(Map.Entry.comparingByKey()).forEach(e -> {
            MolDataValues resData = e.getValue();
            Atom atom = resData.getAtom();
            interpolateRates(atom, field1, field2, field3, 1.0);
        });
    }

    public static void interpolateRates(Atom atom, double field1, double field2, double field3, double tol) {
        var relaxData = atom.getRelaxationData();
        boolean hasDeuterium = relaxData.values().stream().
                anyMatch(v -> v.getRelaxationSet().relaxType() == RelaxTypes.RQ);
        double f1 = Math.pow(field1, -2);
        double f2 = Math.pow(field2, -2);
        double f3 = Math.pow(field3, -2);
        double f = (2 * f2 - f3 - f1) / (2.0 * (f1 - f3));
        if (hasDeuterium) {
            int nTypes = RelaxTypes.values().length;
            Double[][] values = new Double[nTypes][2];
            Double[][] errors = new Double[nTypes][2];
            RelaxationData[][] rData = new RelaxationData[nTypes][2];
            int[] nValues = new int[nTypes];
            for (var entry : relaxData.entrySet()) {
                RelaxationData data = entry.getValue();
                RelaxationSet relaxationSet = entry.getKey();
                double dField = relaxationSet.field();
                int iField = -1;
                if (Math.abs(dField - field1) < tol) {
                    iField = 0;
                } else if (Math.abs(dField - field3) < tol) {
                    iField = 1;
                }
                if (iField != -1) {
                    int iType = relaxationSet.relaxType().ordinal();
                    values[iType][iField] = data.getValue();
                    errors[iType][iField] = data.getError();
                    rData[iType][iField] = data;
                    nValues[iType]++;
                }
            }
            MoleculeBase moleculeBase = MoleculeFactory.getActive();
            var relaxMap = moleculeBase.relaxationSetMap();

            for (int i = 0; i < nTypes; i++) {
                if (nValues[i] == 2 && values[i][0] != null && values[i][1] != null) {
                    double v1 = values[i][0];
                    double v3 = values[i][1];
                    double e1 = errors[i][0];
                    double e3 = errors[i][1];
                    double v2 = v1 * (0.5 + f) + v3 * (0.5 - f);
                    double e2 = e1 * (0.5 + f) + e3 * (0.5 - f);

                    RelaxTypes type = RelaxTypes.values()[i];
                    RelaxationSet relaxationSet1 = rData[i][0].getRelaxationSet();
                    RelaxationSet relaxationSet3 = rData[i][1].getRelaxationSet();
                    ResonanceSource resSource = rData[i][0].getResonanceSource();
                    String id = relaxationSet1.name() + "_interp_" + ((int) field2);
                    RelaxationSet newSet = new RelaxationSet(id, type, field2, relaxationSet1.temperature(), Collections.emptyMap());
                    RelaxationData newData = new RelaxationData(newSet, resSource, v2, e2);
                    resSource.getAtom().addRelaxationData(newSet, newData);
                    relaxationSet1.active(false);
                    relaxationSet3.active(false);
                    relaxMap.put(id, newSet);

                }
            }
        }
    }

}

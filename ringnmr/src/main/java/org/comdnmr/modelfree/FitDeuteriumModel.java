package org.comdnmr.modelfree;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.nmrfx.chemistry.*;
import org.nmrfx.chemistry.relax.RelaxationData;

import java.util.*;

public class FitDeuteriumModel {

    Map<String, MolDataValues> getData(boolean requireCoords) {
        Map<String, MolDataValues> molDataValues = new HashMap<>();
        MoleculeBase mol = MoleculeFactory.getActive();
        System.out.println("mol " + mol);
        if (mol != null) {
            for (Polymer polymer : mol.getPolymers()) {
                System.out.println("poly " + polymer);
                for (Residue residue : polymer.getResidues()) {
                    System.out.println("res " + residue);
                    for (Atom atom : residue.getAtoms()) {
                        var relaxData = atom.getRelaxationData();
                        boolean hasDeuterium = relaxData.values().stream().
                                anyMatch(v -> v.expType == RelaxationData.relaxTypes.RQ);

                        if (hasDeuterium) {
                            Vector3D vec = null;
                            MolDataValues molData = vec == null
                                    ? new MolDataValues(atom)
                                    : new MolDataValues(atom, vec.toArray());
                            Set<Long> fields = new HashSet<>();
                            for (var entry : relaxData.entrySet()) {
                                RelaxationData data = entry.getValue();
                                fields.add(Math.round(data.getField()));
                            }
                            for (var field : fields) {
                                Double r1 = null;
                                Double r1Error = null;

                                Double r2 = null;
                                Double r2Error = null;

                                Double rQ = null;
                                Double rQError = null;

                                Double rAP = null;
                                Double rAPError = null;

                                for (var entry : relaxData.entrySet()) {
                                    RelaxationData data = entry.getValue();
                                    if (Math.round(data.getField()) == field) {
                                        switch (data.getExpType()) {
                                            case R1:
                                                r1 = data.getValue();
                                                r1Error = data.getError();
                                                break;
                                            case R2:
                                                r2 = data.getValue();
                                                r2Error = data.getError();
                                                break;
                                            case RQ:
                                                rQ = data.getValue();
                                                rQError = data.getError();
                                                break;
                                            case RAP:
                                                rAP = data.getValue();
                                                rAPError = data.getError();
                                                break;
                                        }
                                    }
                                }
                                if ((r1 != null) && (r2 != null) && (rQ != null) && (rAP != null)) {
                                    RelaxEquations relaxObj = RelaxEquations.getRelaxEquations(field * 1e6, "H", "N");

                                    RelaxDataValue dValue = new DeuteriumDataValue(molData, r1, r1Error, r2, r2Error,
                                            rQ, rQError, rAP, rAPError, relaxObj);
                                    molData.addData(dValue);
                                }
                            }
                            if (!molData.dataValues.isEmpty()) {
                                molDataValues.put(molData.specifier, molData);
                            }
                        }

                    }
                }
            }
        }
        return molDataValues;
    }

    public void testIsoModel(String searchKey, List<String> modelNames) {
        Map<String, MolDataValues> molData = getData(false);
        if (searchKey != null) {
            if (molData.containsKey(searchKey)) {
                var keepVal = molData.get(searchKey);
                molData.clear();
                molData.put(searchKey, keepVal);
            }
        }

        if (!molData.isEmpty()) {
            testModels(molData, modelNames);
        } else {
            throw new IllegalStateException("No relaxation data to analyze.  Need T1,T2 and NOE");
        }
    }
    public void testModels(Map<String, MolDataValues> molData, List<String> modelNames) {
        Random random = new Random();
        MoleculeBase mol = MoleculeFactory.getActive();
        Map<String, MolDataValues> molDataRes = new TreeMap<>();
        for (String key : molData.keySet()) {
            molDataRes.clear();
            MolDataValues resData = molData.get(key);
            System.out.println(resData);
            if (!resData.getData().isEmpty()) {
                molDataRes.put(key, resData);
                List<Double> rValues = new ArrayList<>();
                List<Double> fields = new ArrayList<>();
                for (var value: resData.dataValues) {
                   var dValue = (DeuteriumDataValue) value;
                    rValues.add(dValue.R1);
                    rValues.add(dValue.R2);
                    rValues.add(dValue.rQ);
                    rValues.add(dValue.rAP);
                    System.out.println("dValue " + dValue.relaxObj.getSF());
                    fields.add(dValue.relaxObj.getSF());
                }
                double[][] mappingResult = DeuteriumMapping.jointMapping(rValues, fields);
                for (int i=0;i<mappingResult[0].length;i++) {
                    double field = mappingResult[0][i];
                    double jValue = mappingResult[1][i];
                    double errValue = mappingResult[2][i];
                    double logJ = Math.log10(jValue);
                    double logJUp = Math.log10(jValue + errValue);
                    double logJLow = Math.log10(jValue - errValue);
                    System.out.printf("%7.1f %7.5f %7.5f %7.5f %7.5f %7.5f\n",field / 1.0e6, jValue, errValue, logJ, logJLow, logJUp);
                }

            }
        }
    }

}

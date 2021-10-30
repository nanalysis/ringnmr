/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.modelfree;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.optim.PointValuePair;
import org.comdnmr.modelfree.RelaxFit.Score;
import org.comdnmr.modelfree.models.MFModelIso;
import org.nmrfx.chemistry.Atom;
import org.nmrfx.chemistry.MoleculeBase;
import org.nmrfx.chemistry.MoleculeFactory;
import org.nmrfx.chemistry.Polymer;
import org.nmrfx.chemistry.relax.RelaxationData;
import org.nmrfx.chemistry.relax.RelaxationRex;
import org.nmrfx.chemistry.Residue;
import org.nmrfx.chemistry.relax.OrderPar;

/**
 *
 * @author brucejohnson
 */
public class FitModel {

    Double tau;
    boolean fitTau = false;
    boolean fitExchange = false;
    double tauFraction = 0.25;
    double lambda = 0.0;
    double t2Limit = 0.0;

    Map<String, MolDataValues> getData(boolean requireCoords) {
        Map<String, MolDataValues> molDataValues = new HashMap<>();
        MoleculeBase mol = MoleculeFactory.getActive();
        if (mol != null) {
            for (Polymer polymer : mol.getPolymers()) {
                for (Residue residue : polymer.getResidues()) {
                    var hAtom = residue.getAtom("H");
                    var nAtom = residue.getAtom("N");
                    if ((hAtom != null) && (nAtom != null)) {
                        var hPoint = hAtom.getPoint();
                        var nPoint = nAtom.getPoint();
                        Vector3D vec = null;
                        if ((hPoint == null) || (nPoint == null)) {
                            if (requireCoords) {
                                continue;
                            }
                        } else {
                            vec = hPoint.subtract(nPoint);
                        }
                        var relaxData = hAtom.getRelaxationData();
                        if ((relaxData == null) || relaxData.isEmpty()) {
                            relaxData = nAtom.getRelaxationData();
                        }
                        if (relaxData != null) {
                            MolDataValues molData = vec == null
                                    ? new MolDataValues(nAtom)
                                    : new MolDataValues(nAtom, vec.toArray());
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

                                Double noe = null;
                                Double noeError = null;

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
                                            case NOE:
                                                noe = data.getValue();
                                                noeError = data.getError();
                                                break;
                                        }
                                    }
                                }
                                if ((r1 != null) && (r2 != null) && (noe != null)) {
                                    String e1 = hAtom.getElementName();
                                    String e2 = nAtom.getElementName();
                                    RelaxEquations relaxObj = RelaxEquations.getRelaxEquations(field * 1e6, "H", "N");

                                    RelaxDataValue dValue = new RelaxDataValue(molData, r1, r1Error, r2, r2Error, noe, noeError, relaxObj);
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

    public void testIsoModel() {
        testIsoModel(null, List.of(1));
    }

    public void testIsoModel(String searchKey, List<Integer> modelNums) {
        Map<String, MolDataValues> molData = getData(false);
        if (searchKey != null) {
            if (molData.containsKey(searchKey)) {
                var keepVal = molData.get(searchKey);
                molData.clear();
                molData.put(searchKey, keepVal);
            }
        }

        if (!molData.isEmpty()) {
            if (tau == null) {
                tau = estimateTau(molData).get("tau");
            }

            testModels(molData, modelNums);
        } else {
            throw new IllegalStateException("No relaxation data to analyze.  Need T1,T2 and NOE");
        }
    }

    public Map<String, Double> estimateTau() {
        Map<String, MolDataValues> molData = getData(false);
        return estimateTau(molData);
    }

    Map<String, Double> estimateTau(Map<String, MolDataValues> molData) {
        Map<Long, List<RelaxDataValue>> map = new HashMap<>();
        for (var entry : molData.entrySet()) {
            MolDataValues value = entry.getValue();
            for (RelaxDataValue rlxValue : value.getData()) {
                long sfMHz = Math.round(rlxValue.relaxObj.getSF() / 1.0e6);
                List<RelaxDataValue> values = map.get(sfMHz);
                if (values == null) {
                    values = new ArrayList<>();
                    map.put(sfMHz, values);
                }
                values.add(rlxValue);
            }
        }
        int max = 0;
        long maxMHz = 0;
        for (long sfMHz : map.keySet()) {
            int size = map.get(sfMHz).size();
            if (size > max) {
                maxMHz = sfMHz;
                max = size;
            }
        }
        if (max > 0) {
            Map<String, Double> tauData = CorrelationTime.estimateTau(maxMHz, "N", map.get(maxMHz));
            return tauData;
        } else {
            return Collections.EMPTY_MAP;
        }
    }

    public void testIsoModel(MFModelIso model, Map<String, MolDataValues> molData, String matchSpec) {
        Map<String, MolDataValues> molDataRes = new TreeMap<>();
        tau = estimateTau(molData).get("tau");
        Map<String, String> extras = new HashMap<>();
        MoleculeBase mol = MoleculeFactory.getActive();

        for (String key : molData.keySet()) {
            if ((matchSpec != null) && !key.equals(matchSpec)) {
                continue;
            }
            molDataRes.clear();
            MolDataValues resData = molData.get(key);
            if (!resData.getData().isEmpty()) {
                resData.setTestModel(model);
                molDataRes.put(key, molData.get(key));
                Score score = tryModel(molDataRes, model, tauFraction, fitTau);
                double[] pars = score.getPars();
                double orderPar = pars[0];
                double rexValue = pars[1];
                RelaxationRex relaxData = new RelaxationRex("order", RelaxationData.relaxTypes.S2,
                        resData.atom, new ArrayList<>(), 700.0, 298, orderPar, 0.0,
                        rexValue, 0.0, extras);
                Atom atom = mol.findAtom(key);
                atom.addRelaxationData("order", relaxData);
            }
        }
    }

    public void testModels(Map<String, MolDataValues> molData, List<Integer> modelNums) {
        Map<String, MolDataValues> molDataRes = new TreeMap<>();
        if (tau == null) {
            tau = estimateTau(molData).get("tau");
        }
        MoleculeBase mol = MoleculeFactory.getActive();

        for (String key : molData.keySet()) {
            molDataRes.clear();
            MolDataValues resData = molData.get(key);

            if (!resData.getData().isEmpty()) {
                molDataRes.put(key, resData);
                MFModelIso bestModel = null;
                Score bestScore = null;
                double lowestAIC = Double.MAX_VALUE;
                for (var modelNum : modelNums) {

                    boolean localFitTau;
                    double localTauFraction;
                    if (overT2Limit(molDataRes, t2Limit)) {
                        localTauFraction = tauFraction;
                        localFitTau = fitTau;
                    } else {
                        localTauFraction = 0.0;
                        localFitTau = false;
                    }
                    System.out.println(localFitTau + " " + localTauFraction);
                    MFModelIso model = MFModelIso.buildModel(modelNum, localFitTau, tau,
                            fitExchange);

                    resData.setTestModel(model);
                    Score score = tryModel(molDataRes, model, localTauFraction, localFitTau);
                    if (score.aic() < lowestAIC) {
                        lowestAIC = score.aic();
                        bestModel = model;
                        bestScore = score;
                    }
                }
                if (bestScore != null) {
                    Atom atom = resData.atom;
                    double[] pars = bestScore.getPars();
                    var parNames = bestModel.getParNames();
                    OrderPar orderPar = new OrderPar(atom, bestScore.rss, bestScore.nValues, "model" + bestModel.getNumber());
                    for (int iPar = 0; iPar < parNames.size(); iPar++) {
                        String parName = parNames.get(iPar);
                        double parValue = pars[iPar];
                        Double parError = null;
                        System.out.println(iPar + " " + parName + " " + parValue);
                        orderPar = orderPar.set(parName, parValue, parError);
                    }
                    if (bestModel.hasTau()) {
                        orderPar = orderPar.set("Tau_e", bestModel.getTau(), 0.0);
                    }

                    orderPar = orderPar.set("model", (double) bestModel.getNumber(), null);
                    if ((bestModel.getNumber() == 5) && (lambda > 1.0e-6)) {
                        orderPar = orderPar.setModel();
                    }
                    System.out.println("order par " + orderPar.toString() + " " + bestScore.rms());
                    atom.addOrderPar("order", orderPar);
                }
            }
        }
    }

    public void testModel(Map<String, MolDataValues> molData, MFModelIso model) {
        Map<String, MolDataValues> molDataRes = new TreeMap<>();
        MoleculeBase mol = MoleculeFactory.getActive();

        for (String key : molData.keySet()) {
            molDataRes.clear();
            MolDataValues resData = molData.get(key);

            if (!resData.getData().isEmpty()) {
                double lowestAIC = Double.MAX_VALUE;
                resData.setTestModel(model);
                molDataRes.put(key, molData.get(key));
                Score score = tryModel(molDataRes, model, tauFraction, fitTau);
                Atom atom = resData.atom;
                double[] pars = score.getPars();
                var parNames = model.getParNames();
                OrderPar orderPar = new OrderPar(atom, score.rss, score.getN(), "model1");
                for (int iPar = 0; iPar < parNames.size(); iPar++) {
                    String parName = parNames.get(iPar);
                    double parValue = pars[iPar];
                    Double parError = null;
                    System.out.println(iPar + " " + parName + " " + parValue);
                    orderPar = orderPar.set(parName, parValue, parError);
                }
                if (!fitTau) {
                    orderPar = orderPar.set("Tau_e", tau, 0.0);
                }
                orderPar = orderPar.setModel();
                System.out.println("order par " + orderPar.toString());
                atom.addOrderPar("order", orderPar);
            }
        }
    }

    boolean overT2Limit(Map<String, MolDataValues> molDataRes, double limit) {
        return molDataRes.values().stream().anyMatch(v -> {
            return v.getData().stream().anyMatch(d -> d.R2 > limit);
        });
    }

    Score tryModel(Map<String, MolDataValues> molDataRes, MFModelIso model, double localTauFraction, boolean localFitTau) {
        RelaxFit relaxFit = new RelaxFit();
        relaxFit.setRelaxData(molDataRes);
        relaxFit.setLambda(lambda);
        model.setTauFraction(localTauFraction);
        double[] start = model.getStart(tau, localFitTau);
        double[] lower = model.getLower(tau, localFitTau);
        double[] upper = model.getUpper(tau, localFitTau);
        PointValuePair fitResult = relaxFit.fitResidueToModel(start, lower, upper);
        var score = relaxFit.score(fitResult.getPoint(), true);
        return score;
    }

    public Double getTau() {
        return tau;
    }

    public void setTau(Double value) {
        tau = value;
    }

    public void setFitTau(boolean value) {
        fitTau = value;
    }

    public void setLambda(double value) {
        this.lambda = value;
    }

    public void setT2Limit(double value) {
        this.t2Limit = value;
    }

    public void setTauFraction(double value) {
        this.tauFraction = value;
    }

}

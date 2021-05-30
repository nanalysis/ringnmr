/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.modelfree;

import java.util.ArrayList;
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
import org.comdnmr.modelfree.models.MFModelIso1;
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

    Map<String, MolDataValues> getData() {
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
                        if ((hPoint != null) && (nPoint != null)) {
                            var relaxData = hAtom.getRelaxationData();
                            if ((relaxData == null) || relaxData.isEmpty()) {
                                relaxData = nAtom.getRelaxationData();
                            }
                            if (relaxData != null) {
                                Vector3D vec = hPoint.subtract(nPoint);
                                MolDataValues molData = new MolDataValues(nAtom, vec.toArray());
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
                                                case T1:
                                                    r1 = data.getValue();
                                                    r1Error = data.getError();
                                                    break;
                                                case T2:
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
        }
        return molDataValues;
    }

    public void testIsoModel() {
        Map<String, MolDataValues> molData = getData();
        if (!molData.isEmpty()) {
            double tau = estimateTau(molData);
            MFModelIso model1e = new MFModelIso1(tau, true);
            MFModelIso model1 = new MFModelIso1(tau, false);
//            MFModelIso model2 = new MFModelIso2(tau, true);
//            MFModelIso model5 = new MFModelIso5(tau, true);
//            MFModelIso model6 = new MFModelIso6(tau, true);
            var models = List.of(model1, model1e);
            testModels(molData, models);
        }
    }

    double estimateTau(Map<String, MolDataValues> molData) {
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
            return tauData.get("tau") * 1.0e-9;
        } else {
            return 5.0e-9;
        }
    }

    public void testIsoModel(MFModelIso model, Map<String, MolDataValues> molData, String matchSpec) {
        Map<String, MolDataValues> molDataRes = new TreeMap<>();
        double tau = estimateTau(molData);
        System.out.println("tau " + tau);
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
                Score score = tryModel(molDataRes, model, tau);
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

    public void testModels(Map<String, MolDataValues> molData, List<MFModelIso> models) {
        Map<String, MolDataValues> molDataRes = new TreeMap<>();
        double tau = estimateTau(molData);
        System.out.println("tau " + tau);
        MoleculeBase mol = MoleculeFactory.getActive();

        for (String key : molData.keySet()) {
            molDataRes.clear();
            MolDataValues resData = molData.get(key);

            if (!resData.getData().isEmpty()) {
                MFModelIso bestModel = null;
                Score bestScore = null;
                double lowestAIC = Double.MAX_VALUE;
                for (var model : models) {
                    resData.setTestModel(model);
                    molDataRes.put(key, molData.get(key));
                    Score score = tryModel(molDataRes, model, tau);
                    if (score.aic() < lowestAIC) {
                        lowestAIC = score.aic();
                        bestModel = model;
                        bestScore = score;
                    }
                }
                if (bestScore != null) {
                    Atom atom = resData.atom;
                    double[] pars = bestScore.getPars();
                    double orderValue = pars[0];
                    Double rexValue = pars.length > 1 ? pars[1] : null;
                    Double rexError = pars.length > 1 ? 0.0 : null;
                    OrderPar orderPar = new OrderPar(atom, orderValue, 0.0, bestScore.rss, "model1");
                    orderPar = orderPar.rEx(rexValue, rexError);
                    atom.addOrderPar("order", orderPar);
                }
            }
        }
    }

    Score tryModel(Map<String, MolDataValues> molDataRes, MFModelIso model, double tau) {
        RelaxFit relaxFit = new RelaxFit();
        relaxFit.setRelaxData(molDataRes);
        double[] start = model.getStart(tau, false);
        double[] lower = model.getLower(tau, false);
        double[] upper = model.getUpper(tau, false);
        PointValuePair fitResult = relaxFit.fitResidueToModel(start, lower, upper);
        var score = relaxFit.score(fitResult.getPoint(), true);
        return score;
    }

}

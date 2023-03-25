/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.modelfree;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Random;
import java.util.Set;
import java.util.TreeMap;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.comdnmr.modelfree.models.MFModelIso;
import org.comdnmr.modelfree.models.MFModelIso2sf;
import org.nmrfx.chemistry.Atom;
import org.nmrfx.chemistry.MoleculeBase;
import org.nmrfx.chemistry.MoleculeFactory;
import org.nmrfx.chemistry.Polymer;
import org.nmrfx.chemistry.relax.*;
import org.nmrfx.chemistry.Residue;

/**
 * @author brucejohnson
 */
public class FitR1R2NOEModel extends FitModel {

    Map<String, MolDataValues> getData(boolean requireCoords) {
        Map<String, MolDataValues> molDataValues = new HashMap<>();
        MoleculeBase mol = MoleculeFactory.getActive();
        if (mol != null) {
            for (Polymer polymer : mol.getPolymers()) {
                for (Residue residue : polymer.getResidues()) {
                    var hAtom = residue.getAtom("H");
                    var nAtom = residue.getAtom("N");
                    if (nAtom == null) {
                        continue;
                    }
                    Vector3D vec = null;
                    if (requireCoords) {
                        if (hAtom != null) {
                            var hPoint = hAtom.getPoint();
                            var nPoint = nAtom.getPoint();
                            if ((hPoint == null) || (nPoint == null)) {
                                continue;
                            } else {
                                vec = hPoint.subtract(nPoint);
                            }
                        } else {
                            continue;
                        }
                    }
                    var relaxData = nAtom.getRelaxationData();

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
                                if (data.getResonanceSource().deleted()) {
                                    break;
                                }
                                if (Math.round(data.getField()) == field) {
                                    switch (data.getExpType()) {
                                        case R1 -> {
                                            r1 = data.getValue();
                                            r1Error = data.getError();
                                        }
                                        case R2 -> {
                                            r2 = data.getValue();
                                            r2Error = data.getError();
                                        }
                                        case NOE -> {
                                            noe = data.getValue();
                                            noeError = data.getError();
                                        }
                                    }
                                }
                            }
                            if ((r1 != null) && (r2 != null) && (noe != null)) {
                                RelaxEquations relaxObj = RelaxEquations.getRelaxEquations(field * 1e6, "H", "N");

                                R1R2NOEDataValue dValue = new R1R2NOEDataValue(molData, r1, r1Error, r2, r2Error, noe, noeError, relaxObj);
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
        return molDataValues;
    }

    public void testIsoModel() {
        testIsoModel(null, List.of("1"));
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
            if (tau == null) {
                tau = estimateTau(molData).get("tau");
            }

            testModels(molData, modelNames);
        } else {
            throw new IllegalStateException("No relaxation data to analyze.  Need T1,T2 and NOE");
        }
    }

    public Optional<OrderPar> testIsoModelResidue(String searchKey, List<String> modelNames) {
        Random random = new Random();
        Map<String, MolDataValues> molData = getData(false);
        MolDataValues molDataValue = molData.get(searchKey);
        return testModels(molDataValue, searchKey, modelNames, random);
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
                List<RelaxDataValue> values = map.computeIfAbsent(sfMHz, k -> new ArrayList<>());
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
            return CorrelationTime.estimateTau(maxMHz, "N", map.get(maxMHz));
        } else {
            return Collections.emptyMap();
        }
    }

    public void testIsoModel(MFModelIso model, Map<String, MolDataValues> molData, String matchSpec) {
        Random random = new Random();
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
                Score score = tryModel(molDataRes, model, tauFraction, fitTau, random);
                double[] pars = score.getPars();
                double orderPar = pars[0];
                double rexValue = pars[1];
                ResonanceSource resSource = new ResonanceSource(resData.atom);
                RelaxationRex relaxData = new RelaxationRex("order", RelaxationData.relaxTypes.S2,
                        resSource, 700.0, 298, orderPar, 0.0,
                        rexValue, 0.0, extras);
                Atom atom = mol.findAtom(key);
                atom.addRelaxationData("order", relaxData);
            }
        }
    }

    public void testModels(Map<String, MolDataValues> molData, List<String> modelNames) {
        Random random = new Random();
        if (tau == null) {
            tau = estimateTau(molData).get("tau");
        }
        MoleculeBase mol = MoleculeFactory.getActive();
        Map<String, MolDataValues> molDataRes = new TreeMap<>();
        for (String key : molData.keySet()) {
            molDataRes.clear();
            MolDataValues resData = molData.get(key);
            if (!resData.getData().isEmpty()) {
                molDataRes.put(key, resData);
                testModels(resData, key, modelNames, random);
            }
        }
    }

    public Optional<OrderPar> testModels(MolDataValues resData, String key, List<String> modelNames, Random random) {
        Optional<OrderPar> result = Optional.empty();
        Map<String, MolDataValues> molDataRes = new TreeMap<>();
        molDataRes.put(key, resData);

        MFModelIso bestModel = null;
        Score bestScore = null;
        double lowestAIC = Double.MAX_VALUE;
        boolean localFitTau;
        double localTauFraction;
        if (overT2Limit(molDataRes, t2Limit)) {
            localTauFraction = tauFraction;
            localFitTau = fitTau;
        } else {
            localTauFraction = 0.0;
            localFitTau = false;
        }
        for (var modelName : modelNames) {

            MFModelIso model = MFModelIso.buildModel(modelName,
                    localFitTau, tau, localTauFraction, fitExchange);

            resData.setTestModel(model);
            Score score = tryModel(molDataRes, model, localTauFraction, localFitTau, random);
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
            double[][] repData = null;
            if (nReplicates > 2) {
                repData = replicates(molDataRes, bestModel, localTauFraction, localFitTau, pars, random);
            }
            ResonanceSource resSource = new ResonanceSource(resData.atom);

            OrderPar orderPar = new OrderPar(resSource, bestScore.rss, bestScore.nValues, bestScore.nPars, bestModel.getName());
            for (int iPar = 0; iPar < parNames.size(); iPar++) {
                String parName = parNames.get(iPar);
                double parValue = pars[iPar];
                Double parError = null;
                if (repData != null) {
                    DescriptiveStatistics sumStat = new DescriptiveStatistics(repData[iPar]);
                    parError = sumStat.getStandardDeviation();
                }
                orderPar = orderPar.set(parName, parValue, parError);
            }
            if (!bestModel.fitTau()) {
                orderPar = orderPar.set("Tau_e", bestModel.getTau(), 0.0);
            }

            orderPar = orderPar.set("model", (double) bestModel.getNumber(), null);
            if ((bestModel instanceof MFModelIso2sf) && (lambda > 1.0e-6)) {
                orderPar = orderPar.setModel();
            }
            atom.addOrderPar("order", orderPar);
            double[][] jValues = resData.getJValues();
            SpectralDensity spectralDensity = new SpectralDensity(key, jValues);
            atom.addSpectralDensity(key, spectralDensity);
            result = Optional.of(orderPar);
        }
        return result;
    }

    public void testModel(Map<String, MolDataValues> molData, MFModelIso model) {
        Map<String, MolDataValues> molDataRes = new TreeMap<>();
        Random random = new Random();

        for (String key : molData.keySet()) {
            molDataRes.clear();
            MolDataValues resData = molData.get(key);

            if (!resData.getData().isEmpty()) {
                resData.setTestModel(model);
                molDataRes.put(key, molData.get(key));
                Score score = tryModel(molDataRes, model, tauFraction, fitTau, random);
                Atom atom = resData.atom;
                double[] pars = score.getPars();
                var parNames = model.getParNames();
                ResonanceSource resSource = new ResonanceSource(resData.atom);

                OrderPar orderPar = new OrderPar(resSource, score.rss, score.getN(), score.nPars, "model1");
                for (int iPar = 0; iPar < parNames.size(); iPar++) {
                    String parName = parNames.get(iPar);
                    double parValue = pars[iPar];
                    Double parError = null;
                    orderPar = orderPar.set(parName, parValue, parError);
                }
                if (!fitTau) {
                    orderPar = orderPar.set("Tau_e", tau, 0.0);
                }
                orderPar = orderPar.setModel();
                atom.addOrderPar("order", orderPar);
            }
        }
    }

    boolean overT2Limit(Map<String, MolDataValues> molDataRes, double limit) {
        return molDataRes.values().stream().anyMatch(v -> v.getData().stream().anyMatch(d -> d.R2 > limit));
    }

    Score tryModel(Map<String, MolDataValues> molDataRes, MFModelIso model, double localTauFraction, boolean localFitTau, Random random) {
        RelaxFit relaxFit = new RelaxFit();
        relaxFit.setRelaxData(molDataRes);
        relaxFit.setLambda(lambda);
        relaxFit.setFitJ(fitJ);
        model.setTauFraction(localTauFraction);
        double[] start = model.getStart();
        double[] lower = model.getLower();
        double[] upper = model.getUpper();
        double[] keepStart = start.clone();
        int nTries = 3;
        PointValuePair best = null;
        for (int i = 0; i < nTries; i++) {
            PointValuePair fitResult = relaxFit.fitResidueToModel(start, lower, upper);
            System.out.println(i + " testModel " + fitResult.getValue());
            if ((i == 0) || (fitResult.getValue() < best.getValue())) {
                best = fitResult;
            }
            for (int j = 0; j < start.length; j++) {
                start[j] = keepStart[j] + random.nextGaussian() * 0.1 * (upper[j] - lower[j]);
            }
        }
        var score = relaxFit.score(best.getPoint(), true);
        System.out.println(score.rms() + " " + score.value(lambda) + " " + score.complexity() + " " + score.parsOK());
        for (var d : score.getPars()) {
            System.out.print(d + " ");
        }
        System.out.println();
        return score;
    }

    public void setTauFraction(double value) {
        this.tauFraction = value;
    }

    private void appendValueError(StringBuilder stringBuilder, Double val, Double err, String format) {
        if (val != null) {
            stringBuilder.append(String.format(format, val)).append("\t");
        } else {
            stringBuilder.append("\t\t");
        }
        if (err != null) {
            stringBuilder.append(String.format(format, err)).append("\t");
        } else {
            stringBuilder.append("\t\t");
        }
    }

    public void writeData(File file) throws IOException {
        try (FileWriter fileWriter = new FileWriter(file)) {
            var molDataValues = getData(false);
            String header = "residue\tfield\tr1\tr1_err\tr2\tr2_err\tnoe\tnoe_err\n";
            fileWriter.write(header);
            for (var molData : molDataValues.values()) {
                var data = molData.getData();
                for (var value : data) {
                    R1R2NOEDataValue dValue = (R1R2NOEDataValue) value;
                    StringBuilder stringBuilder = new StringBuilder();
                    stringBuilder.append(molData.specifier).append("\t");
                    stringBuilder.append(String.format("%.2f", dValue.relaxObj.getSF() / 1.0e6)).append("\t");
                    appendValueError(stringBuilder, dValue.R1, dValue.R1err, "%.3f");
                    appendValueError(stringBuilder, dValue.R2, dValue.R2err, "%.3f");
                    appendValueError(stringBuilder, dValue.NOE, dValue.NOEerr, "%.3f");
                    fileWriter.write(stringBuilder.toString());
                    fileWriter.write("\n");
                }
            }
        }
    }
}

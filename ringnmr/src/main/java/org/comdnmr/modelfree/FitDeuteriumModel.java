package org.comdnmr.modelfree;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.comdnmr.modelfree.models.MFModelIso;
import org.comdnmr.modelfree.models.MFModelIso2sf;
import org.nmrfx.chemistry.*;
import org.nmrfx.chemistry.relax.OrderPar;
import org.nmrfx.chemistry.relax.RelaxationData;
import org.nmrfx.chemistry.relax.ResonanceSource;
import org.nmrfx.chemistry.relax.SpectralDensity;

import java.util.*;

public class FitDeuteriumModel extends FitModel {

    public static Map<String, MolDataValues> getData(boolean requireCoords) {
        Map<String, MolDataValues> molDataValues = new HashMap<>();
        MoleculeBase mol = MoleculeFactory.getActive();
        if (mol != null) {
            for (Polymer polymer : mol.getPolymers()) {
                for (Residue residue : polymer.getResidues()) {
                    for (Atom atom : residue.getAtoms()) {
                        MolDataValues molData = getMolDataValues(atom);
                        if ((molData != null) && !molData.dataValues.isEmpty()) {
                            molDataValues.put(molData.specifier, molData);
                        }
                    }
                }
            }
        }
        return molDataValues;
    }

    public static boolean hasDeuteriumData(Atom atom) {
        var relaxData = atom.getRelaxationData();
        boolean hasDeuterium = relaxData.values().stream().
                anyMatch(v -> v.expType == RelaxationData.relaxTypes.RQ);
        return hasDeuterium;
    }

    public static MolDataValues getMolDataValues(Atom atom) {
        var relaxData = atom.getRelaxationData();
        boolean hasDeuterium = relaxData.values().stream().
                anyMatch(v -> v.expType == RelaxationData.relaxTypes.RQ);
        MolDataValues molData = null;
        if (hasDeuterium) {
            Vector3D vec = null;
            molData = vec == null
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
                    RelaxEquations relaxObj = RelaxEquations.getRelaxEquations(field * 1e6, "D", "C");

                    RelaxDataValue dValue = new DeuteriumDataValue(molData, r1, r1Error, r2, r2Error,
                            rQ, rQError, rAP, rAPError, relaxObj);
                    molData.addData(dValue);
                }
            }
        }
        return molData;
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

    public void calcJ(Map<String, MolDataValues> molData) {
        Map<String, MolDataValues> molDataRes = new TreeMap<>();
        for (String key : molData.keySet()) {
            molDataRes.clear();
            MolDataValues resData = molData.get(key);
            if (!resData.getData().isEmpty()) {
                molDataRes.put(key, resData);
                List<Double> rValues = new ArrayList<>();
                List<Double> errValues = new ArrayList<>();
                List<Double> fields = new ArrayList<>();
                for (var value : resData.dataValues) {
                    var dValue = (DeuteriumDataValue) value;
                    rValues.add(dValue.R1);
                    rValues.add(dValue.R2);
                    rValues.add(dValue.rQ);
                    rValues.add(dValue.rAP);
                    errValues.add(dValue.R1err);
                    errValues.add(dValue.R2err);
                    errValues.add(dValue.rQError);
                    errValues.add(dValue.rAPError);
                    fields.add(dValue.relaxObj.getSF());
                }
                double[][] mappingResult = DeuteriumMapping.jointMapping(rValues, errValues, fields);
                for (int i = 0; i < mappingResult[0].length; i++) {
                    double field = mappingResult[0][i];
                    double jValue = mappingResult[1][i];
                    double errValue = mappingResult[2][i];
                    double logJ = Math.log10(jValue);
                    double logJUp = Math.log10(jValue + errValue);
                    double logJLow = Math.log10(jValue - errValue);
                }

            }
        }
    }

    public void testModels(Map<String, MolDataValues> molData, List<String> modelNames) {
        Random random = new Random();
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
        boolean localFitTau = fitTau;
        double localTauFraction = tauFraction;
        for (var modelName : modelNames) {

            MFModelIso model = MFModelIso.buildModel(modelName,
                    localFitTau, tau, localTauFraction, fitExchange);

            resData.setTestModel(model);
            Score score = tryModel(molDataRes, model, localTauFraction, localFitTau, random);
            if ((score != null) && (score.aic() < lowestAIC)) {
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
            atom.addOrderPar("order_parameter_list_1", orderPar);
            double[][] jValues  = resData.getJValues();
            SpectralDensity spectralDensity = new SpectralDensity(key,jValues);
            atom.addSpectralDensity(key, spectralDensity);

            result = Optional.of(orderPar);
        }
        return result;
    }

    Score tryModel(Map<String, MolDataValues> molDataRes, MFModelIso model, double localTauFraction, boolean localFitTau, Random random) {
        RelaxFit relaxFit = new RelaxFit();
        relaxFit.setRelaxData(molDataRes);
        relaxFit.setLambda(0.0);
        relaxFit.setFitJ(true);
        model.setTauFraction(localTauFraction);
        double[] start = model.getStart();
        double[] lower = model.getLower();
        double[] upper = model.getUpper();
        double[] keepStart = start.clone();
        var parNames = model.getParNames();
        for (var parName : parNames) {
            System.out.println(parName);
        }
        for (int i = 0; i < start.length; i++) {
            System.out.println(i + " start " + start[i] + " lower " + lower[i] + " upper " + upper[i]);
        }
        int nTries = 3;
        PointValuePair best = null;
        for (int i = 0; i < nTries; i++) {
            try {
                PointValuePair fitResult = relaxFit.fitResidueToModel(start, lower, upper);
                if ((i == 0) || (fitResult.getValue() < best.getValue())) {
                    best = fitResult;
                }
            } catch(Exception iaE) {
                System.out.println(iaE.getMessage());
            }
            for (int j = 0; j < start.length; j++) {
                start[j] = keepStart[j] + random.nextGaussian() * 0.1 * (upper[j] - lower[j]);
            }
        }
        if (best != null) {
            var score = relaxFit.score(best.getPoint(), true);
            System.out.println(score.rms() + " " + score.value() + " " + score.complexity() + " " + score.parsOK());
            System.out.print("Pars result: ");
            for (var d : score.getPars()) {
                System.out.print(d + " ");
            }
            System.out.println();
            return score;
        }
        return null;
    }

    double[][] replicates(Map<String, MolDataValues> molDataRes,
                          MFModelIso bestModel, double localTauFraction,
                          boolean localFitTau, double[] pars, Random random) {
        double[][] repData = new double[pars.length][nReplicates];
        for (int iRep = 0; iRep < nReplicates; iRep++) {
            boolean ok = false;
            for (int jTry = 0; jTry < 5; jTry++) {
                Score score2 = fitReplicate(molDataRes, bestModel, localTauFraction, localFitTau, pars, random);
                if (score2 != null) {
                    double[] repPars = score2.getPars();
                    for (int iPar = 0; iPar < pars.length; iPar++) {
                        repData[iPar][iRep] = repPars[iPar];
                    }
                    ok = true;
                    break;
                }
            }
            if (!ok) {
                return null;
            }
        }
        return repData;
    }

    Score fitReplicate(Map<String, MolDataValues> molDataRes, MFModelIso model,
                       double localTauFraction, boolean localFitTau, double[] pars, Random random) {
        RelaxFit relaxFit = new RelaxFit();
        relaxFit.setRelaxData(molDataRes);
        relaxFit.setLambda(lambda);
        relaxFit.setFitJ(fitJ);
        Map<String, MolDataValues> molDataMap = relaxFit.genBootstrap(random, model, pars);
        relaxFit.setRelaxData(molDataMap);

        model.setTauFraction(localTauFraction);
        double[] lower = model.getLower();
        double[] upper = model.getUpper();
        PointValuePair fitResult = relaxFit.fitResidueToModel(pars, lower, upper);
        if (fitResult == null) {
            return null;
        } else {
            return relaxFit.score(fitResult.getPoint(), true);
        }
    }

}

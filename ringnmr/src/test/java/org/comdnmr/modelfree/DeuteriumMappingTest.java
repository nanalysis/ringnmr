package org.comdnmr.modelfree;

import org.comdnmr.data.DynamicsSource;
import org.comdnmr.modelfree.models.MFModelIso1;
import org.comdnmr.modelfree.models.MFModelIso1f;
import org.comdnmr.modelfree.models.MFModelIso2f;
import org.junit.Assert;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class DeuteriumMappingTest {
    List<Double> fields = List.of(4.0E8, 5.0E8, 8.0E8, 9.0E8);

    List<Double> rValuesWithErrs = List.of(
            26.36059572, 0.404816211, 106.1659187, 1.819686576, 24.61990335, 1.331138025, 88.85133494, 3.474549981,
            22.33003784, 0.148052025, 105.3207282, 2.126010315, 19.87279621, 0.368971996, 94.31636064, 1.047352748,
            16.76840591, 0.207212294, 101.560506, 1.471771902, 13.32686669, 0.216933607, 90.1311929, 1.364259199,
            15.95096119, 0.236173794, 97.48258269, 0.666038741, 12.37328619, 0.311354286, 90.2230769, 1.396329853);

    @Test
    public void testDeuteriumFit1Field() {
        List<Double> rValues = new ArrayList<>();
        List<Double> rValueErrs = new ArrayList<>();
        for (int i = 0; i < rValuesWithErrs.size(); i += 2) {
            rValues.add(rValuesWithErrs.get(i));
            rValueErrs.add(rValuesWithErrs.get(i + 1));
        }
        for (int j = 0; j < fields.size(); j++) {
            var rvaluesField = rValues.subList(j * 4, j * 4 + 4);
            var rvalueErrsField = rValueErrs.subList(j * 4, j * 4 + 4);
            RelaxEquations rlxEq = RelaxEquations.getRelaxEquations(fields.get(j), "D", "C");
            var field400 = List.of(rlxEq.getW()[1]);
            var jValues = DeuteriumMapping.jointMapping(rvaluesField, rvalueErrsField, field400);
            double r1 = rlxEq.R1_D(jValues[1]);
            double r2 = rlxEq.R2_D(jValues[1]);
            double rQ = rlxEq.RQ_D(jValues[1]);
            double rAP = rlxEq.Rap_D(jValues[1]);
            double[] result = {r1, r2, rQ, rAP};
            for (int i = 0; i < result.length; i++) {
                result[i] = (result[i] - rvaluesField.get(i)) / rValueErrs.get(i);
            }
            double[] zeros = {0.0, 0.0, 0.0, 0.0};
            Assert.assertArrayEquals(zeros, result, 1.0);
        }
    }

    @Test
    public void testDeuteriumFit() {
        List<Double> rValues = new ArrayList<>();
        List<Double> rValueErrs = new ArrayList<>();
        for (int i = 0; i < rValuesWithErrs.size(); i += 2) {
            rValues.add(rValuesWithErrs.get(i));
            rValueErrs.add(rValuesWithErrs.get(i + 1));
        }
        var jValues = DeuteriumMapping.jointMapping(rValues, rValueErrs, fields);
        double[] jField = new double[3];
        for (int j = 0; j < fields.size(); j++) {
            RelaxEquations rlxEq = RelaxEquations.getRelaxEquations(fields.get(j), "D", "C");
            System.arraycopy(jValues[4], j * 3, jField, 0, 3);
            double r1 = rlxEq.R1_D(jField);
            double r2 = rlxEq.R2_D(jField);
            double rQ = rlxEq.RQ_D(jField);
            double rAP = rlxEq.Rap_D(jField);
            double[] result = {r1, r2, rQ, rAP};
            var rvaluesField = rValues.subList(j * 4, j * 4 + 4);
            var rvalueErrsField = rValueErrs.subList(j * 4, j * 4 + 4);
            for (int i = 0; i < result.length; i++) {
                result[i] = (result[i] - rvaluesField.get(i)) / rValueErrs.get(i);
            }
            double[] zeros = {0.0, 0.0, 0.0, 0.0};
            Assert.assertArrayEquals(zeros, result, 3.0);
        }
    }

    @Test
    public void testModel1f() {
        RelaxEquations rlxEq = RelaxEquations.getRelaxEquations(5.0e8, "D", "C");
        var model = new MFModelIso1f();
        model.setSScale(9.0);
        double[] pars = {6.2, 0.7, 0.04};
        var jValues = model.calc(rlxEq.wValues, pars);
        double r1 = rlxEq.R1_D(jValues);
        double r2 = rlxEq.R2_D(jValues);
        double rQ = rlxEq.RQ_D(jValues);
        double rAP = rlxEq.Rap_D(jValues);
        List<Double> rValues = List.of(r1, r2, rQ, rAP);
        List<Double> rValuesErrs = List.of(r1*0.04, r2*0.03, rQ*0.05, rAP*0.02);

        var fields = List.of(rlxEq.getW()[1]);
        var jValuesCalc = DeuteriumMapping.jointMapping(rValues, rValuesErrs, fields);
        Assert.assertArrayEquals(jValues, jValuesCalc[1], 1.0e-12);
    }

    @Test
    public void testModel2f() {
        RelaxEquations rlxEq = RelaxEquations.getRelaxEquations(9.0e8, "D", "C");
        var model = new MFModelIso2f();
        model.setSScale(9.0);
        double[] pars = {8.2, 0.7, 0.04, 1.0};
        var jValues = model.calc(rlxEq.wValues, pars);
        double r1 = rlxEq.R1_D(jValues);
        double r2 = rlxEq.R2_D(jValues);
        double rQ = rlxEq.RQ_D(jValues);
        double rAP = rlxEq.Rap_D(jValues);
        List<Double> rValues = List.of(r1, r2, rQ, rAP);
        List<Double> rValuesErrs = List.of(r1*0.04, r2*0.03, rQ*0.05, rAP*0.02);

        var fields = List.of(rlxEq.getW()[1]);
        var jValuesCalc = DeuteriumMapping.jointMapping(rValues, rValuesErrs, fields);
        Assert.assertArrayEquals(jValues, jValuesCalc[1], 1.0e-12);
    }
    @Test
    public void testModel1fFit() {
        DynamicsSource dynamicsSourceFactory = new DynamicsSource(true, true, true, true);

        List<Double> rValues = new ArrayList<>();
        List<Double> rValueErrs = new ArrayList<>();
        double[] v = {0.0, 1.0, 2.0};
        MolDataValues resData = new MolDataValues("3.CB",v, dynamicsSourceFactory);
        for (int i = 0; i < rValuesWithErrs.size(); i += 8) {
            double r1 = rValuesWithErrs.get(i);
            double r1Error = rValuesWithErrs.get(i+1);
            double r2 = rValuesWithErrs.get(i+2);
            double r2Error = rValuesWithErrs.get(i+3);
            double rQ = rValuesWithErrs.get(i+4);
            double rQError = rValuesWithErrs.get(i+5);
            double rAP = rValuesWithErrs.get(i+6);
            double rAPError = rValuesWithErrs.get(i+7);
            RelaxEquations relaxObj = RelaxEquations.getRelaxEquations(fields.get(i/8), "D", "C");


            rValues.add(rValuesWithErrs.get(i));
            rValueErrs.add(rValuesWithErrs.get(i + 1));
            RelaxDataValue dValue = new DeuteriumDataValue(resData, r1, r1Error, r2, r2Error,
                    rQ, rQError, rAP, rAPError, relaxObj);
            resData.addData(dValue);
        }
        FitDeuteriumModel fitModel = new FitDeuteriumModel();
        fitModel.setTau(12.0);
        fitModel.setFitTau(true);
        fitModel.setTauFraction(0.5);
        fitModel.setNReplicates(10);
        fitModel.setFitJ(true);
        Random random = new Random();
        var modelNames = List.of("D1f");
        var result = fitModel.testModels(resData, "tst", modelNames, random);
        Assert.assertTrue(result.isPresent());
        Assert.assertEquals(0.0, result.get().orderPar().getReducedChiSqr(), 3.0);
    }

}
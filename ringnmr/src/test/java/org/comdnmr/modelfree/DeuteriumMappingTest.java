package org.comdnmr.modelfree;

import org.comdnmr.modelfree.models.MFModelIso1;
import org.comdnmr.modelfree.models.MFModelIso1f;
import org.comdnmr.modelfree.models.MFModelIso2f;
import org.junit.Assert;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

public class DeuteriumMappingTest {
    @Test
    public void testIndependentMapping() {
        double R1 = 3.0;
        double R1rho = 9.0;
        double RQ = 3.0;
        double Rap = 1.0;
        var jValues = DeuteriumMapping.independentMapping(R1, R1rho, RQ, Rap);
        for (var jValue : jValues) {
            System.out.print(jValue + " ");
        }
        System.out.println();
        Assert.assertEquals(1, 1, 1);
    }

    @Test
    public void testIndependentMapping2() {

        List<Double> rValuesWithErrs = List.of(
                26.36059572, 0.404816211, 106.1659187, 1.819686576, 24.61990335, 1.331138025, 88.85133494, 3.474549981,
                22.33003784, 0.148052025, 105.3207282, 2.126010315, 19.87279621, 0.368971996, 94.31636064, 1.047352748,
                16.76840591, 0.207212294, 101.560506, 1.471771902, 13.32686669, 0.216933607, 90.1311929, 1.364259199,
                15.95096119, 0.236173794, 97.48258269, 0.666038741, 12.37328619, 0.311354286, 90.2230769, 1.396329853);

        List<Double> rValues = new ArrayList<>();
        List<Double> rValueErrs = new ArrayList<>();
        for (int i = 0; i < rValuesWithErrs.size(); i += 2) {
            rValues.add(rValuesWithErrs.get(i));
            rValueErrs.add(rValuesWithErrs.get(i + 1));
        }

        List<Double> fields = List.of(4.0E8, 5.0E8, 8.0E8, 9.0E8);


        var jValues = DeuteriumMapping.jointMapping(rValues, rValueErrs, fields);
        for (var jValue : jValues) {
            System.out.print(jValue + " ");
        }
        System.out.println();
        Assert.assertEquals(1, 1, 1);
    }

    @Test
    public void testDeuteriumFit() {
        List<Double> fields = List.of(4.0E8, 5.0E8, 8.0E8, 9.0E8);

        List<Double> rValuesWithErrs = List.of(
                26.36059572, 0.404816211, 106.1659187, 1.819686576, 24.61990335, 1.331138025, 88.85133494, 3.474549981,
                22.33003784, 0.148052025, 105.3207282, 2.126010315, 19.87279621, 0.368971996, 94.31636064, 1.047352748,
                16.76840591, 0.207212294, 101.560506, 1.471771902, 13.32686669, 0.216933607, 90.1311929, 1.364259199,
                15.95096119, 0.236173794, 97.48258269, 0.666038741, 12.37328619, 0.311354286, 90.2230769, 1.396329853);

        List<Double> rValues = new ArrayList<>();
        List<Double> rValueErrs = new ArrayList<>();
        for (int i = 0; i < rValuesWithErrs.size(); i += 2) {
            rValues.add(rValuesWithErrs.get(i));
            rValueErrs.add(rValuesWithErrs.get(i + 1));
        }
        var rvalues400 = rValues.subList(0,4);
        var rvalueErrs400 = rValueErrs.subList(0,4);
        RelaxEquations rlxEq = RelaxEquations.getRelaxEquations(4.0e8,"D", "C");
        var field400 = List.of( rlxEq.getW()[1]);
        var jValues = DeuteriumMapping.jointMapping(rvalues400,rvalueErrs400, field400);
        for (int i=0;i<jValues.length;i++) {
            for (int j=0;j<jValues[i].length;j++) {
                System.out.print(jValues[i][j] + " " );
            }
            System.out.println();
        }
        double r1 = rlxEq.R1_D(jValues[1]);
        double r2 = rlxEq.R2_D(jValues[1]);
        double rQ = rlxEq.RQ_D(jValues[1]);
        double rAP = rlxEq.Rap_D(jValues[1]);
        System.out.println("r1 " + r1 + " r2 " + r2 + " rQ " + rQ + " rAP " + rAP);
    }

    @Test
    public void testModel1f() {
        RelaxEquations rlxEq = RelaxEquations.getRelaxEquations(9.0e8,"D", "C");
        var model = new MFModelIso1f();
        model.setSScale(9.0);
        double[] pars = {8.2,0.7, 0.04};
        var jValues = model.calc(rlxEq.wValues,pars);
        for (var jValue:jValues) {
            System.out.println(jValue);
        }
        double r1 = rlxEq.R1_D(jValues);
        double r2 = rlxEq.R2_D(jValues);
        double rQ = rlxEq.RQ_D(jValues);
        double rAP = rlxEq.Rap_D(jValues);
        System.out.println("r1 " + r1 + " r2 " + r2 + " rQ " + rQ + " rAP " + rAP);

    }
    @Test
    public void testModel2f() {
        RelaxEquations rlxEq = RelaxEquations.getRelaxEquations(9.0e8,"D", "C");
        var model = new MFModelIso2f();
        model.setSScale(9.0);
        double[] pars = {8.2,0.7, 0.04, 1.0};
        var jValues = model.calc(rlxEq.wValues,pars);
        for (var jValue:jValues) {
            System.out.println(jValue);
        }
        double r1 = rlxEq.R1_D(jValues);
        double r2 = rlxEq.R2_D(jValues);
        double rQ = rlxEq.RQ_D(jValues);
        double rAP = rlxEq.Rap_D(jValues);
        System.out.println("r1 " + r1 + " r2 " + r2 + " rQ " + rQ + " rAP " + rAP);

    }
}
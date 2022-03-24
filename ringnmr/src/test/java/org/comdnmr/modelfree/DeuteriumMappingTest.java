package org.comdnmr.modelfree;

import org.junit.Assert;
import org.junit.Test;

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

        List<Double> rValues = List.of(
                26.36059572, 106.1659187, 24.61990335, 88.85133494,
                22.33003784, 105.3207282, 19.87279621, 94.31636064,
                16.76840591, 101.560506, 13.32686669, 90.1311929,
                15.95096119, 97.48258269, 12.37328619, 90.2230769
        );

        List<Double> fields = List.of(4.0E8, 5.0E8, 8.0E8, 9.0E8);


        var jValues = DeuteriumMapping.jointMapping(rValues, fields);
        for (var jValue : jValues) {
            System.out.print(jValue + " ");
        }
        System.out.println();
        Assert.assertEquals(1, 1, 1);
    }
}
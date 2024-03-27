package org.comdnmr.modelfree;

import org.junit.Assert;
import org.junit.Test;

public class RelaxEquationsTest {
    double twoPI = 2.0 * Math.PI;
    double[] R1rhoTargets = {
        367.95510639514123, 296.2386769001938, 227.1116787400174, 171.8855136375552, 131.6300869486979,
        71.688811251761, 71.09001993233547, 70.37120833078939, 69.53968236184136, 68.60363518294635,
        7.334645374625591, 7.333994404029964, 7.333198930185943, 7.332259047151564, 7.331174866020584
    };
    double[] R1rhoTaus = {1.0e-5, 1.0e-6, 1.0e-7};
    double[] R1rhoOmegaEs = {twoPI * 20.0e3, twoPI * 25.0e3, twoPI * 30.0e3, twoPI * 35.0e3, twoPI * 40.0e3};
    double R1rhoOmegaR = twoPI * 10.0e3;
    double R1rhoS2 = 0.3295;
    double R1rhoSIGMA = -37.0e-6;

    @Test
    public void testR1rhoCSA() {
        RelaxEquations.setSigma("C", R1rhoSIGMA);
        RelaxEquations relax = new RelaxEquations(400.0e6, "H", "C");
        int i = 0;
        for (double R1rhoTau: R1rhoTaus) {
            for (double R1rhoOmegaE: R1rhoOmegaEs) {
               double R1rhoTarget = R1rhoTargets[i];
                Assert.assertEquals(
                    R1rhoTarget,
                    relax.R1rhoCSA(
                       R1rhoOmegaR, R1rhoOmegaE, R1rhoTau, R1rhoS2
                    ),
                    1.0e-6
                );
                i++;
            }
        }
    }
}
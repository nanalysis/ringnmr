/*
 * CoMD/NMR Software : A Program for Analyzing NMR Dynamics Data
 * Copyright (C) 2018-2019 Bruce A Johnson
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
 /*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.modelfree.models;

import java.util.List;

/**
 *
 * @author brucejohnson
 */
public class MFModelIso5 extends MFModelIso4 {

    double tauS;
    double complexity = 0.0;

    public MFModelIso5() {
        super();
        nPars = 4;
    }

    public MFModelIso5(double tauM) {
        super(tauM);
        nPars = 4;
    }

    public MFModelIso5(boolean includeEx) {
        super(includeEx);
        nPars = includeEx ? 5 : 4;
    }

    public MFModelIso5(double tauM, boolean includeEx) {
        super(tauM, includeEx);
        nPars = includeEx ? 5 : 4;
    }

    @Override
    public List<String> getParNames() {
        return getAllParNames("S2", "Tau_f", "Sf2", "Tau_s");
    }

    @Override
    public double[] calc(double[] omegas) {
        double tauMx = tauM * 1.0e-9;
        double tauFx = tauF * 1.0e-9;
        double tauSx = tauS * 1.0e-9;
        double[] J = new double[omegas.length];
        int j = 0;
        double ss2 = s2 / sf2;
        complexity = 0.0;
        for (double omega : omegas) {
            double omega2 = omega * omega;
            double vM = s2 * tauMx / (1.0 + omega2 * tauMx * tauMx);
            double vMS = sf2 * (1.0 - ss2) * (tauMx * tauSx * (tauMx + tauSx))
                    / (tauMx * tauMx * tauSx * tauSx * omega2 + (tauMx + tauSx) * (tauMx + tauSx));
            double vMF = (1.0 - sf2) * ss2 * (tauMx * tauFx * (tauMx + tauFx))
                    / (tauMx * tauMx * tauFx * tauFx * omega2 + (tauMx + tauFx) * (tauMx + tauFx));
            double vMFS = (1.0 - sf2) * (1.0 - ss2) * (tauFx * tauMx * tauSx * (tauFx * (tauMx + tauSx) + tauMx * tauSx))
                    / (tauFx * tauFx * tauMx * tauMx * tauSx * tauSx * omega2
                    + (tauFx * (tauMx + tauSx) + tauMx * tauSx) * (tauFx * (tauMx + tauSx) + tauMx * tauSx));
            // complexity += Math.abs(vMS / vM) + Math.abs(vMF / vM) + Math.abs(vMFS / vM);
            complexity
                    += Math.abs(1.0 - sf2)
                    + Math.abs(1.0 - ss2)
                    + 0.1 * Math.abs(1.0 - s2)
                    + tauSx * 1.0e9
                    + tauFx * 1.0e9;
            if (hasTau) {
                //   complexity += Math.abs(7.75 + Math.log10(tauM))
            }
            J[j++] = 0.4 * (vM + vMS + vMF + vMFS);
        }
        return J;
    }

    @Override
    public double[] calc(double[] omegas, double[] pars) {
        int parStart = 0;
        if (!hasTau) {
            tauM = pars[0];
            parStart = 1;
        }

        this.s2 = pars[parStart];
        this.tauF = pars[parStart + 1];
        this.sf2 = pars[parStart + 2];
        this.tauS = pars[parStart + 3];
        return calc(omegas);
    }

    @Override
    public double getComplexity() {
        return complexity;
    }

    public double[] calc(double[] omegas, double s2, double tauF, double sf2, double tauS) {
        this.s2 = s2;
        this.tauF = tauF;
        this.sf2 = sf2;
        this.tauS = tauS;
        return calc(omegas);
    }

    @Override
    public boolean checkParConstraints() {
        return tauS * 1.0e-9 < tauM && tauF < tauS;
    }

    @Override
    public double calcWeight() {
        return 0.0;
    }

    @Override
    public double[] getStart(double tau, boolean includeTau) {
        if (includeEx) {
            return getParValues(includeTau, tau, 0.9, 0.1, 0.9, 0.01, 2.0);
        } else {
            return getParValues(includeTau, tau, 0.9, 0.1, 0.9, 0.01);
        }
    }

    @Override
    public double[] getLower(double tau, boolean includeTau) {
        if (includeEx) {
            return getParValues(includeTau, tauLower(tau), 0.0, 0.001, 0.0, 0.15, 0.0);
        } else {
            return getParValues(includeTau, tauLower(tau), 0.0, 0.001, 0.0, 0.15);
        }
    }

    @Override
    public double[] getUpper(double tau, boolean includeTau) {
        if (includeEx) {
            return getParValues(includeTau, tauUpper(tau), 1.0, 0.15, 1.0, tau / 2.0, 100.0);
        } else {
            return getParValues(includeTau, tauUpper(tau), 1.0, 0.15, 1.0, tau / 2.0);

        }
    }

    @Override
    public int getNumber() {
        return 5;
    }

}

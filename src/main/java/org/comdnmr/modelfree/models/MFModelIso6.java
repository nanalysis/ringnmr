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
public class MFModelIso6 extends MFModelIso5 {

    double tauS;
    double complexity = 0.0;

    public MFModelIso6() {
        super();
        nPars = 4;
    }

    public MFModelIso6(double tauM) {
        super(tauM);
        nPars = 4;
    }

    public MFModelIso6(boolean includeEx) {
        super(includeEx);
        nPars = includeEx ? 5 : 4;
    }

    public MFModelIso6(double tauM, boolean includeEx) {
        super(tauM, includeEx);
        nPars = includeEx ? 5 : 4;
    }

    @Override
    public List<String> getParNames() {
        return getAllParNames("S2", "Tau_f", "Sf2", "Tau_s");
    }

    @Override
    public double[] calc(double[] omegas) {
        double[] J = new double[omegas.length];
        int j = 0;
        double ss2 = s2 / sf2;
        complexity = 0.0;
        for (double omega : omegas) {
            double omega2 = omega * omega;
            double vM = s2 * tauM / (1.0 + omega2 * tauM * tauM);
            double vMS = sf2 * (1.0 - ss2) * (tauM * tauS * (tauM + tauS))
                    / (tauM * tauM * tauS * tauS * omega2 + (tauM + tauS) * (tauM + tauS));
            double vMF = (1.0 - sf2) * ss2 * (tauM * tauF * (tauM + tauF))
                    / (tauM * tauM * tauF * tauF * omega2 + (tauM + tauF) * (tauM + tauF));
            double vMFS = (1.0 - sf2) * (1.0 - ss2) * (tauF * tauM * tauS * (tauF * (tauM + tauS) + tauM * tauS))
                    / (tauF * tauF * tauM * tauM * tauS * tauS * omega2
                    + (tauF * (tauM + tauS) + tauM * tauS) * (tauF * (tauM + tauS) + tauM * tauS));
            // complexity += Math.abs(vMS / vM) + Math.abs(vMF / vM) + Math.abs(vMFS / vM);
            complexity
                    += Math.abs(1.0 - sf2)
                    + Math.abs(1.0 - ss2)
                    + 0.1 * Math.abs(1.0 - s2)
                    + Math.abs(7.75 + Math.log10(tauM))
                    + Math.abs(9.825 + Math.log10(tauS))
                    + Math.abs(12.0 + Math.log10(tauF));
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
        this.tauF = pars[parStart + 1] * 1.0e-9;
        this.sf2 = pars[parStart + 2];
        this.tauS = pars[parStart + 3] * 1.0e-9;
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
            return getParValues(includeTau, tau, 0.9, 0.05, 0.9, 0.3, 2.0);
        } else {
            return getParValues(includeTau, tau, 0.9, 0.05, 0.9, 0.5);
        }
    }

    @Override
    public double[] getLower(double tau, boolean includeTau) {
        if (includeEx) {
            return getParValues(includeTau, tau / 5., 0.0, 1.0e-3, 0.0, 1.0e-3, 0.0);
        } else {
            return getParValues(includeTau, tau / 5., 0.0, 1.0e-3, 0.0, 0.15);
        }
    }

    @Override
    public double[] getUpper(double tau, boolean includeTau) {
        if (includeEx) {
            return getParValues(includeTau, tau * 5., 1.0, 0.15, 1.0, tau * 1.0e9 / 2.0, 100.0);
        } else {
            return getParValues(includeTau, tau * 5., 1.0, 0.15, 1.0, tau * 1.0e9 / 2.0);

        }
    }

    @Override
    public int getNumber() {
        return 6;
    }

}

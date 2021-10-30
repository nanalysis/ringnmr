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
public class MFModelIso2s extends MFModelIso1s {

    double sf2;

    public MFModelIso2s() {
        super();
        nPars = 3;
    }

    public MFModelIso2s(double tauM) {
        super(tauM);
        nPars = 3;
    }

    public MFModelIso2s(boolean includeEx) {
        super(includeEx);
        nPars = includeEx ? 4 : 3;
    }

    public MFModelIso2s(double tauM, boolean includeEx) {
        super(tauM, includeEx);
        nPars = includeEx ? 4 : 3;
    }

    @Override
    public List<String> getParNames() {
        return getAllParNames("Sf2", "Tau_s", "Ss2");
    }

    @Override
    public double[] calc(double[] omegas) {
        double tauMx = tauM * 1.0e-9;
        double tauSx = tauS * 1.0e-9;
        double[] J = new double[omegas.length];
        int j = 0;
        double s2 = sf2 * ss2;
        for (double omega : omegas) {
            double omega2 = omega * omega;
            double taus = tauMx * tauSx / (tauMx + tauSx);
            double value1 = s2 * tauMx / (1.0 + omega2 * tauMx * tauMx);
            double value2 = (sf2 - s2) * (taus) / (1.0 + omega2 * taus * taus);
            J[j++] = 0.4 * (value1 + value2);
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

        this.sf2 = pars[parStart];
        this.tauS = pars[parStart + 1];
        this.ss2 = pars[parStart + 2];
        return calc(omegas);
    }

    public double[] calc(double[] omegas, double s2, double tauS, double sf2) {
        this.sf2 = s2;
        this.tauS = tauS;
        this.ss2 = sf2;
        return calc(omegas);
    }

    @Override
    public boolean checkParConstraints() {
        return tauS < tauM;
    }

    @Override
    public double[] getStart(double tau, boolean includeTau) {
        if (includeEx) {
            return getParValues(includeTau, tau, 0.9, tau / 5.0, 0.9, 2.0);
        } else {
            return getParValues(includeTau, tau, 0.9, tau / 5.0, 0.9);
        }
    }

    @Override
    public double[] getLower(double tau, boolean includeTau) {
        if (includeEx) {
            return getParValues(includeTau, tauLower(tau), 0.0, SLOW_LIMIT, 0.0, 0.0);
        } else {
            return getParValues(includeTau, tauLower(tau), 0.0, SLOW_LIMIT, 0.0);
        }
    }

    @Override
    public double[] getUpper(double tau, boolean includeTau) {
        if (includeEx) {
            return getParValues(includeTau, tauUpper(tau), 1.0, tau / 2.0, 1.0, 100.0);
        } else {
            return getParValues(includeTau, tauUpper(tau), 1.0, tau / 2.0, 1.0);

        }
    }

    @Override
    public int getNumber() {
        return 5;
    }

}

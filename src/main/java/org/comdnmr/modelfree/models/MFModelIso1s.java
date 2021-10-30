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
public class MFModelIso1s extends MFModelIso {

    double ss2;
    double tauS;

    public MFModelIso1s() {
        super(false);
        nPars = 2;
    }

    public MFModelIso1s(double tauM) {
        super(true);
        this.tauM = tauM;
        nPars = 2;
    }

    public MFModelIso1s(boolean includeEx) {
        super(false, includeEx);
        nPars = includeEx ? 3 : 2;
    }

    public MFModelIso1s(double tauM, boolean includeEx) {
        super(true, includeEx);
        this.tauM = tauM;
        nPars = includeEx ? 3 : 2;
    }

    @Override
    public List<String> getParNames() {
        return getAllParNames("Ss2", "Tau_s");
    }

    @Override
    public double[] calc(double[] omegas) {
        double tauMx = tauM * 1.0e-9;
        double tauSx = tauS * 1.0e-9;
        double[] J = new double[omegas.length];
        int j = 0;
        for (double omega : omegas) {
            double omega2 = omega * omega;
            double tauf = tauMx * tauSx / (tauMx + tauSx);
            double value1 = ss2 * tauMx / (1.0 + omega2 * tauMx * tauMx);
            double value2 = (1.0 - ss2) * (tauf) / (1.0 + omega2 * tauf * tauf);
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

        this.ss2 = pars[parStart];
        this.tauS = pars[parStart + 1];
        return calc(omegas);
    }

    public double[] calc(double[] omegas, double s2, double tauS) {
        this.ss2 = s2;
        this.tauS = tauS;
        return calc(omegas);
    }

    @Override
    public boolean checkParConstraints() {
        return tauS < tauM;
    }

    @Override
    public double[] getStart(double tau, boolean includeTau) {
        if (includeEx) {
            return getParValues(includeTau, tau, 0.9, tau / 5.0, 2.0);
        } else {
            return getParValues(includeTau, tau, 0.9, tau / 5.0);
        }
    }

    @Override
    public double[] getLower(double tau, boolean includeTau) {
        if (includeEx) {
            return getParValues(includeTau, tauLower(tau), 0.0, SLOW_LIMIT, 0.0);
        } else {
            return getParValues(includeTau, tauLower(tau), 0.0, SLOW_LIMIT);
        }
    }

    @Override
    public double[] getUpper(double tau, boolean includeTau) {
        if (includeEx) {
            return getParValues(includeTau, tauUpper(tau), 1.0, tau / 2.0, 100.0);
        } else {
            return getParValues(includeTau, tauUpper(tau), 1.0, tau / 2.0);

        }
    }

    @Override
    public int getNumber() {
        return 3;
    }

}

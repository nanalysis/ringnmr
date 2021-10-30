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

import java.util.ArrayList;
import java.util.List;
import org.comdnmr.modelfree.RelaxFit;

/**
 *
 * @author brucejohnson
 */
public class MFModelIso1 extends MFModelIso {

    double s2;

    public MFModelIso1() {
        super(false);
        nPars = 1;
    }

    public MFModelIso1(boolean includeEx) {
        super(false, includeEx);
        nPars = includeEx ? 2 : 1;
    }

    public MFModelIso1(double tauM) {
        super(true);
        this.tauM = tauM;
        nPars = 1;
    }

    public MFModelIso1(double tauM, boolean includeEx) {
        super(true, includeEx);
        this.tauM = tauM;
        nPars = includeEx ? 2 : 1;
    }

    @Override
    public List<String> getParNames() {
        return getAllParNames("S2");
    }

    @Override
    public double[] calc(double[] omegas) {
        double tauMx  = tauM * 1.0e-9;
        double[] J = new double[omegas.length];
        int j = 0;
        for (double omega : omegas) {
            double omega2 = omega * omega;
            J[j++] = 0.4 * s2 * tauMx / (1.0 + omega2 * tauMx * tauMx);
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
        return calc(omegas);
    }

    public double[] calc(double[] omegas, double s2) {
        this.s2 = s2;
        return calc(omegas);
    }

    @Override
    public double[] getStart(double tau, boolean includeTau) {
        if (includeEx) {
            return getParValues(includeTau, tau, 0.9, 2.0);
        } else {
            return getParValues(includeTau, tau, 0.9);
        }
    }

    @Override
    public double[] getLower(double tau, boolean includeTau) {
        if (includeEx) {
            return getParValues(includeTau, tauLower(tau), 0.0, 0.0);
        } else {
            return getParValues(includeTau, tauLower(tau), 0.0);
        }
    }

    @Override
    public double[] getUpper(double tau, boolean includeTau) {
        if (includeEx) {
            return getParValues(includeTau, tauUpper(tau), 1.0, 100.0);
        } else {
            return getParValues(includeTau, tauUpper(tau), 1.0);

        }
    }

    @Override
    public int getNumber() {
        return 1;
    }

}

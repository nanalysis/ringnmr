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
import org.comdnmr.modelfree.RelaxFit;

/**
 *
 * @author brucejohnson
 */
public class MFModelIso4 extends MFModelIso2 {

    double sf2;

    public MFModelIso4() {
        super();
        nPars = 3;
    }

    public MFModelIso4(double tauM) {
        super(tauM);
        nPars = 3;
    }

    public MFModelIso4(boolean includeEx) {
        super(includeEx);
        nPars = includeEx ? 4 : 3;
    }

    public MFModelIso4(double tauM, boolean includeEx) {
        super(tauM, includeEx);
        nPars = includeEx ? 4 : 3;
    }

    @Override
    public List<String> getParNames() {
        return getAllParNames("S2", "Tau_f", "Sf2");
    }

    @Override
    public double[] calc(double[] omegas) {
        double tauMx = tauM * 1.0e-9;
        double tauFx = tauF * 1.0e-9;
        double[] J = new double[omegas.length];
        int j = 0;
        for (double omega : omegas) {
            double omega2 = omega * omega;
            double tauf = tauMx * tauFx / (tauMx + tauFx);
            double value1 = s2 * tauMx / (1.0 + omega2 * tauMx * tauMx);
            double value2 = (sf2 - s2) * (tauf) / (1.0 + omega2 * tauf * tauf);
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

        this.s2 = pars[parStart];
        this.tauF = pars[parStart + 1];
        this.sf2 = pars[parStart + 2];
        return calc(omegas);
    }

    public double[] calc(double[] omegas, double s2, double tauF, double sf2) {
        this.s2 = s2;
        this.tauF = tauF;
        this.sf2 = sf2;
        return calc(omegas);
    }

    @Override
    public boolean checkParConstraints() {
        return tauF < tauM;
    }

    @Override
    public double[] getStart(double tau, boolean includeTau) {
        if (includeEx) {
            return getParValues(includeTau, tau, 0.9, tau / 40.0, 0.9, 2.0);
        } else {
            return getParValues(includeTau, tau, 0.9, tau / 40.0, 0.9);
        }
    }

    @Override
    public double[] getLower(double tau, boolean includeTau) {
        if (includeEx) {
            return getParValues(includeTau, tauLower(tau), 0.0, tau / 1000.0, 0.0, 0.0);
        } else {
            return getParValues(includeTau, tauLower(tau), 0.0, tau / 1000.0, 0.0);
        }
    }

    @Override
    public double[] getUpper(double tau, boolean includeTau) {
        if (includeEx) {
            return getParValues(includeTau, tauUpper(tau), 1.0, tau / 10.0, 1.0, 100.0);
        } else {
            return getParValues(includeTau, tauUpper(tau), 1.0, tau / 10.0, 1.0);

        }
    }

    @Override
    public int getNumber() {
        return 4;
    }

}

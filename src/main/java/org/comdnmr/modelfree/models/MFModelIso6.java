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
        for (double omega : omegas) {
            double omega2 = omega * omega;
            double taue = tauM * tauS / (tauM + tauS);
            double tauf = tauM * tauF / (tauM + tauF);
            double value1 = s2 * tauM / (1.0 + omega2 * tauM * tauM);
            double value2 = (sf2 - s2) * (taue) / (1.0 + omega2 * taue * taue);//(eS[i]*eS[i] + w2);
            double value3 = (1.0 - sf2) * (tauf) / (1.0 + omega2 * tauf * tauf);//(eF[i]*eF[i] + w2);
            J[j++] = 0.4 * (value1 + value2 + value3);
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

    public double[] calc(double[] omegas, double s2, double tauF, double sf2, double tauS) {
        this.s2 = s2;
        this.tauF = tauF;
        this.sf2 = sf2;
        this.tauS = tauS;
        return calc(omegas);
    }

    @Override
    public boolean checkParConstraints() {
        return tauF < tauM && tauS < tauM;
    }

    @Override
    public double calcWeight() {
        return 0.0;
    }

    @Override
    public double[] getStart(double tau, boolean includeTau) {
        if (includeEx) {
            return getParValues(includeTau, tau, 0.9, tau / 40.0, 0.9, tau / 40.0, 2.0);
        } else {
            return getParValues(includeTau, tau, 0.9, tau / 40.0, 0.9, tau / 40.0);
        }
    }

    @Override
    public double[] getLower(double tau, boolean includeTau) {
        if (includeEx) {
            return getParValues(includeTau, tau / 10., 0.0, tau / 1000.0, 0.0, tau / 1000.0, 0.0);
        } else {
            return getParValues(includeTau, tau / 10., 0.0, tau / 1000.0, 0.0, tau / 1000.0);
        }
    }

    @Override
    public double[] getUpper(double tau, boolean includeTau) {
        if (includeEx) {
            return getParValues(includeTau, tau * 10., 1.0, tau / 10.0, 1.0, tau / 10.0, 100.0);
        } else {
            return getParValues(includeTau, tau * 10., 1.0, tau / 10.0, 1.0, tau / 10.0);

        }
    }

}

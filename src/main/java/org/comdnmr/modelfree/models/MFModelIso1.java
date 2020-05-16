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

import org.comdnmr.modelfree.RelaxFit;

/**
 *
 * @author brucejohnson
 */
public class MFModelIso1 extends MFModelIso {

    double tauM;
    double s2;
    final boolean hasTau;

    public MFModelIso1() {
        hasTau = false;
    }

    public MFModelIso1(double tauM) {
        hasTau = true;
        this.tauM = tauM;
    }

    @Override
    public double[] calc(double[] omegas) {
        double[] J = new double[omegas.length];
        int j = 0;
        for (double omega : omegas) {
            double omega2 = omega * omega;
            J[j++] = 0.4 * s2 * tauM / (1.0 + omega2 * tauM * tauM);
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

}

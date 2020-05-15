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
public class MFModel6 extends MFModel5 {

    double tauS;

    public MFModel6(RelaxFit.DiffusionType diffType, double[][] D, double[][] VT, double[] v) {
        super(diffType, D, VT, v);
    }

    @Override
    double calc(double omega2, int i) {
        double[] eF = diffPars.calcDiffusione(tauF);
        double[] eS = diffPars.calcDiffusione(tauS);
        double value1 = s2 * (Df[i] * a[i]);
        double value2 = (sf2 - s2) * (eS[i] * a[i]) / (1.0 + omega2 * eS[i] * eS[i]);
        double value3 = (1.0 - sf2) * (eF[i] * a[i]) / (1.0 + omega2 * eF[i] * eF[i]);
        return value1 + value2 + value3;

    }

    @Override
    public double[] calc(double[] omegas, double[] pars) {
        this.s2 = pars[0];
        this.tauF = pars[1];
        this.sf2 = pars[2];
        this.tauS = pars[3];
        return calc(omegas);
    }

    public double[] calc(double[] omegas, double s2, double tauF, double sf2, double tauS) {
        this.s2 = s2;
        this.tauF = tauF;
        this.sf2 = sf2;
        this.tauS = tauS;
        return calc(omegas);
    }

}

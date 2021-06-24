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
import org.comdnmr.modelfree.DiffusionPars;
import org.comdnmr.modelfree.RelaxFit;

/**
 *
 * @author brucejohnson
 */
public abstract class MFModelAniso extends MFModel {

    DiffusionPars diffPars;
    double[] dDiff;
    double[] a;
    double[] Df;

    MFModelAniso(RelaxFit.DiffusionType diffType,
            double[] v) {
        diffPars = new DiffusionPars(diffType, v);
    }

    MFModelAniso(RelaxFit.DiffusionType diffType,
            double[][] D, double[][] VT, double[] v) {
        diffPars = new DiffusionPars(diffType, D, VT, v);
        dDiff = diffPars.dDiff;
        a = diffPars.a;
    }

    public List<String> getAllParNames(String... pars) {
        var parNames = new ArrayList<String>();
        for (var par : pars) {
            parNames.add(par);
        }
        if (includeEx) {
            parNames.add("Rex");
        }
        return parNames;
    }

    public void update(double[][] D, double[][] VT) {
        diffPars.update(D, VT);
        dDiff = diffPars.dDiff;
        a = diffPars.a;
    }

    abstract double calc(double omega2, int i);

    public double[] calc(double[] omegas) {
        double[] J = new double[omegas.length];
        int j = 0;
        for (double omega : omegas) {
            double omega2 = omega * omega;
            Df = diffPars.getDf(omega2);
            double sum = 0.0;
            for (int i = 0; i < dDiff.length; i++) {
                sum += calc(omega2, i);
            }
            J[j++] = 0.4 * sum;
        }
        return J;
    }

}

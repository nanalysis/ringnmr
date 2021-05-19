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

/**
 *
 * @author brucejohnson
 */
public abstract class MFModel {

    int nPars;
    boolean includeEx;

    public abstract double[] calc(double[] omegas, double[] pars);

    public abstract double[] calc(double[] omega);

    public double[] getParValues(boolean includeTau, double... parValues) {
        int n = includeTau ? nPars + 1 : nPars;
        int start = includeTau ? 0 : 1;
        double[] values = new double[n];
        System.arraycopy(parValues, start, values, 0, n);
        return values;

    }

    public abstract double[] getStart(double tau, boolean includeTau);

    public abstract double[] getLower(double tau, boolean includeTau);

    public abstract double[] getUpper(double tau, boolean includeTau);

    public boolean checkParConstraints() {
        return true;
    }

    public boolean includesEx() {
        return includeEx;
    }

    public int getNPars() {
        return nPars;
    }

}

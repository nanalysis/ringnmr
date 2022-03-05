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

/**
 *
 * @author brucejohnson
 */
public abstract class MFModelIso extends MFModel {

    double tauM;
    double rEX;
    double tauFrac = 0.25;

    public MFModelIso(boolean fitTau, double targetTau, double tauFraction,
            boolean includeEx) {
        this.fitTau = fitTau;
        this.targetTau = targetTau;
        this.tauFrac = tauFrac;
        this.includeEx = includeEx;
        if (!fitTau) {
            tauM = targetTau;
        }
    }

    public MFModelIso(double targetTau) {
        this(false, targetTau, 0.0, false);
    }

    public MFModelIso() {
        this(true, 0.0, 0.0, false);
    }

    public List<String> getAllParNames(String... pars) {
        var parNames = new ArrayList<String>();
        if (fitTau) {
            parNames.add("Tau_e");
        }
        for (var par : pars) {
            parNames.add(par);
        }
        if (includeEx) {
            parNames.add("Rex");
        }
        return parNames;
    }

    public boolean fitTau() {
        return fitTau;
    }

    public double getTau() {
        return tauM;
    }

    public void setTauFraction(double value) {
        tauFrac = value;
    }

    public double tauLower() {
        return targetTau - targetTau * tauFrac;
    }

    public double tauUpper() {
        return targetTau + targetTau * tauFrac;
    }

    public abstract double[] getStart();

    public static MFModelIso buildModel(String modelName, boolean fitTau,
            double tau, double tauFrac,
            boolean fitExchange) {
        MFModelIso model;
        switch (modelName) {
            case "1":
                model = new MFModelIso1(fitTau, tau, tauFrac, fitExchange);
                break;
            case "1f":
                model = new MFModelIso1f(fitTau, tau, tauFrac, fitExchange);
                break;
            case "1s":
                model = new MFModelIso1s(fitTau, tau, tauFrac, fitExchange);
                break;
            case "2s":
                model = new MFModelIso2s(fitTau, tau, tauFrac, fitExchange);
                break;
            case "2f":
                model = new MFModelIso2f(fitTau, tau, tauFrac, fitExchange);
                break;
            case "2sf":
                model = new MFModelIso2sf(fitTau, tau, tauFrac, fitExchange);
                break;
            default:
                model = null;
        }
        return model;

    }

}

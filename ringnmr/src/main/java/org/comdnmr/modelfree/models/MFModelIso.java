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
import java.util.Arrays;
import java.util.List;

import static org.comdnmr.modelfree.models.MFModelIso2sf.TAU_PRIME;

/**
 * @author brucejohnson
 */
public abstract class MFModelIso extends MFModel {

    private static final String[] MODEL_NAMES = {"1", "1f", "1s", "2s", "2sf", "1sf"};

    double sN = 1.0;
    double tauM;
    double rEX;
    double tauFrac;

    public MFModelIso(boolean fitTau, double targetTau, double tauFraction,
                      boolean includeEx) {
        this.fitTau = fitTau;
        this.targetTau = targetTau;
        this.tauFrac = tauFraction;
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

    public static String[] getAllModelNames() {
        return MODEL_NAMES;
    }

    public List<String> getAllParNames(String... pars) {
        var parNames = new ArrayList<String>();
        if (fitTau) {
            parNames.add("Tau_e");
        }
        parNames.addAll(Arrays.asList(pars));
        if (includeEx) {
            parNames.add("Rex");
        }
        return parNames;
    }

    public void setSScale(double value) {
        sN = value;
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
    /** Enforce a canonical labelling: a lone internal mode lives in the fast
     *  channel, and when both are active tauF < tauS. The likelihood is invariant
     *  under this swap, so it changes only the labels. */
    static void canonicalise(MFModelIso2sf model) {
        double sf2 = model.getSf2();
        double tauF = model.getTauF();
        double ss2 = model.getSs2();
        double tauS = model.getTauS();
        double tauM = model.getTau();

        boolean fastOff = (tauF <= 0.0) || (1.0 - sf2 / model.getSN() <= 1e-6);
        boolean slowOn  = (tauS >  0.0) && (1.0 - ss2 > 1e-6);

        boolean swap = (fastOff && slowOn)                   // lone mode -> fast
                || (!fastOff && slowOn && tauF > tauS);      // both on -> order them
        double[] pars;
        int start;
        if (model.fitTau) {
            pars = new double[5];
            pars[0] = tauM;
            start = 1;
        } else {
            pars = new double[4];
            start = 0;
        }
        if (swap) {
            pars[start] = ss2;
            pars[start + 1] = tauS;
            pars[start + 2] = sf2;
            pars[start + 3] = tauF;
        } else {
            pars[start] = sf2;
            pars[start + 1] = tauF;
            pars[start + 2] = ss2;
            pars[start + 3] = tauS;
        }
    }

    protected double[] createStandardPars(double sf2, double tauF, double ss2, double tauS) {
        double[] pars;

        double sN = 1.0;
        if (this instanceof MFModelIso2sf mfModelIso2sf) {
            sN = mfModelIso2sf.getSN();
        }

        boolean fastOff = (tauF <= 0.0) || (1.0 - sf2 / sN <= 1e-6);
        boolean slowOn  = (tauS >  0.0) && (1.0 - ss2 > 1e-6);

        boolean swap = (fastOff && slowOn)                   // lone mode -> fast
                || (!fastOff && slowOn && tauF > tauS);      // both on -> order them
        int start;
        if (fitTau) {
            pars = new double[5];
            pars[0] = tauM;
            start = 1;
        } else {
            pars = new double[4];
            start = 0;
        }
        if (swap) {
            pars[start] = ss2;
            pars[start + 1] = tauS;
            pars[start + 2] = sf2;
            pars[start + 3] = tauF;
        } else {
            pars[start] = sf2;
            pars[start + 1] = tauF;
            pars[start + 2] = ss2;
            pars[start + 3] = tauS;
        }

        return pars;
    }

    public abstract double[] getStandardPars(double[] pars);

    public static MFModelIso buildModel(String modelName, boolean fitTau,
                                        double tau, double tauFrac,
                                        boolean fitExchange) {
        MFModelIso model;
        if (modelName.startsWith("model")) {
            modelName = modelName.substring(5);
        }
        model = switch (modelName) {
            case "1", "D1" -> new MFModelIso1(fitTau, tau, tauFrac, fitExchange);
            case "1f", "D1f" -> new MFModelIso1f(fitTau, tau, tauFrac, fitExchange);
            case "1s", "D1s" -> new MFModelIso1s(fitTau, tau, tauFrac, fitExchange);
            case "1sf", "D1sf" -> new MFModelIso1sf(fitTau, tau, tauFrac, fitExchange);
            case "2s", "D2s" -> new MFModelIso2s(fitTau, tau, tauFrac, fitExchange);
            case "2sf", "D2sf" -> new MFModelIso2sf(fitTau, tau, tauFrac, fitExchange);
            case "2sfx", "D2sfx" -> new MFModelIso2sf(fitTau, tau, tauFrac, fitExchange);
            default -> throw new IllegalArgumentException("Unknown model " + modelName);
        };
        if (modelName.charAt(0) == 'D') {
            model.setSScale(9.0);
        }
        return model;
    }
}

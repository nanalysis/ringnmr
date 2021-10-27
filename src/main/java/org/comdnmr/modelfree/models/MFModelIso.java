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
    final boolean hasTau;

    public MFModelIso(boolean hasTau) {
        this(hasTau, false);
    }

    public MFModelIso(boolean hasTau, boolean includeEx) {
        this.hasTau = hasTau;
        this.includeEx = includeEx;
    }

    public List<String> getAllParNames(String... pars) {
        var parNames = new ArrayList<String>();
        if (!hasTau) {
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

    public static MFModelIso buildModel(int modelNum, boolean fitTau, double tau,
            boolean fitExchange) {
        MFModelIso model;
        switch (modelNum) {
            case 1:
                model = fitTau ? new MFModelIso1(fitExchange)
                        : new MFModelIso1(tau, fitExchange);
                break;
            case 2:
                model = fitTau ? new MFModelIso2(fitExchange)
                        : new MFModelIso2(tau, fitExchange);
                break;
            case 5:
                model = fitTau ? new MFModelIso5(fitExchange)
                        : new MFModelIso5(tau, fitExchange);
                break;
            case 6:
                model = fitTau ? new MFModelIso6(fitExchange)
                        : new MFModelIso6(tau, fitExchange);
                break;
            default:
                model = null;
        }
        return model;

    }

}

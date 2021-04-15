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
package org.comdnmr.data;

import java.util.HashMap;

/**
 *
 * @author Bruce Johnson
 */
public class CPMGExperiment extends Experiment {
    
    double[] fieldStrengths;

    public CPMGExperiment(String name, String nucleus, double field, double temperature) {
        super(name, nucleus, field, temperature, "CPMG");
        this.fieldStrengths = fieldStrengths;
    }

    public void setXVals(double[] xVals) {
        this.fieldStrengths = xVals;
    }

    public double[] getXVals() {
        return fieldStrengths;
    }

}

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
package org.comdnmr.eqnfit;

/**
 *
 * @author Bruce Johnson
 */
public class ParValue implements ParValueInterface {

    String state;
    String residue;
    String resName;
    String name;
    double value;
    double err;

    public ParValue(String parName) {
        this.name = parName;
    }

    public ParValue(String parName, double value) {
        this.name = parName;
        this.value = value;
    }

    public ParValue(String residue, String state, String parName, double value, double err) {
        this.name = parName;
        this.value = value;
        this.err = err;
    }

    @Override
    public String getName() {
        return name;
    }

    @Override
    public double getValue() {
        return value;
    }

    @Override
    public double getError() {
        return err;
    }

    @Override
    public String getResidue() {
        return residue;
    }
    
    @Override
    public String getResidueName() {
        return resName;
    }

    @Override
    public String getState() {
        return state;
    }

}

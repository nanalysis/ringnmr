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
package org.comdnmr.gui;

/**
 * Data extra values for storing value, lowPercentile and highPercentile.
 */
public class ErrorExtraValues {

    private double value;
    private double lowPercentile;
    private double highPercentile;

    public ErrorExtraValues(double value, double lowPercentile, double highPercentile) {
        this.value = value;
        this.lowPercentile = lowPercentile;
        this.highPercentile = highPercentile;
    }

    public double getValue() {
        return value;
    }

    public double getLowPercentile() {
        return lowPercentile;
    }

    public double getHighPercentile() {
        return highPercentile;
    }

    @Override
    public String toString() {
        return value + " " + lowPercentile + " " + highPercentile;
    }
}

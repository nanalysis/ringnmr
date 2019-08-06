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
 * @author teddycolon
 */
public class CESTPeak {

    double position;
    double depth;
    double width50;
    double width50LB;
    double width50UB;
    double width25;
    double width25LB;
    double width25UB;
    double width75;
    double width75LB;
    double width75UB;
    double baseline;
    int pkInd;

    CESTPeak(int peakInd, double... peakInfo) {
        this.pkInd = peakInd;
        this.position = peakInfo[0];
        this.depth = peakInfo[1];
        this.width50 = peakInfo[2];
        this.width50LB = peakInfo[3];
        this.width50UB = peakInfo[4];
        this.width25 = peakInfo[5];
        this.width25LB = peakInfo[6];
        this.width25UB = peakInfo[7];
        this.width75 = peakInfo[8];
        this.width75LB = peakInfo[9];
        this.width75UB = peakInfo[10];
        this.baseline = peakInfo[11];
    }

    public double[] getHalfWidths() {
        double[] halfWidths = new double[3];
        halfWidths[0] = (width25LB > width25UB) ? width25UB : width25LB;
        halfWidths[1] = (width50LB > width50UB) ? width50UB : width50LB;
        halfWidths[2] = (width75LB > width75UB) ? width75UB : width75LB;
        
        return halfWidths;
    }
            
    public double getDepth() {
        return depth;
    }

    public double getPosition() {
        return position;
    }

    public Double[] getWidths() {
        Double[] widths = new Double[3];
        widths[0] = width25; // 25 %
        widths[1] = width50; // 50 %
        widths[2] = width75; // 75 %
        return widths;
    }

}

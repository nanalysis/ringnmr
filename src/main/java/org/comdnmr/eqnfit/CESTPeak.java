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

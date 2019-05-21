/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.fit.calc;

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

//    CESTPeak(double position, double depth, double width, double widthLB, double widthUB, double width2, double widthLB2, double widthUB2, double width3, double widthLB3, double widthUB3, int pkInd) {
//        this.position = position;
//        this.depth = depth;
//        this.width50 = width;
//        this.width50LB = widthLB;
//        this.width50UB = widthUB;
//        this.width25 = width2;
//        this.width25LB = widthLB2;
//        this.width25UB = widthUB2;
//        this.width75 = width3;
//        this.width75LB = widthLB3;
//        this.width75UB = widthUB3;
//        this.pkInd = pkInd;
//    }
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

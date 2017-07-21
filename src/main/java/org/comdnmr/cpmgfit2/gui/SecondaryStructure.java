/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.cpmgfit2.gui;

import javafx.scene.paint.Color;

/**
 *
 * @author brucejohnson
 */
public class SecondaryStructure {

    int start;
    int end;
    int mode;
    String label;
    Color color = Color.GRAY;

    public SecondaryStructure(int start, int end, int mode, String label) {
        this.start = start;
        this.end = end;
        this.mode = mode;
        this.label = label;
        if (this.label.equals(".")) {
            this.label = "";
        }
    }

    public SecondaryStructure(int start, int end, String strMode, String label, String colorName) {
        this.start = start;
        this.end = end;
        if (strMode.startsWith("h") || strMode.startsWith("H")) {
            this.mode = 1;
        } else if (strMode.startsWith("s") || strMode.startsWith("S")) {
            this.mode = 2;
        } else if (strMode.startsWith("d") || strMode.startsWith("D")) {
            this.mode = 3;
        } else {
            this.mode = 0;
        }
        this.label = label;
        if (this.label.equals(".")) {
            this.label = "";
        }
        this.color = Color.web(colorName);
    }

    boolean isHelix() {
        return mode == 1;
    }

    boolean isSheet() {
        return mode == 2;
    }

    boolean isCoil() {
        return mode == 0;
    }

    boolean isDisulfide() {
        return mode == 3;
    }

    int getStart() {
        return start;
    }

    int getEnd() {
        return end;
    }

    Color getColor() {
        return color;
    }

    String getLabel() {
        return label;
    }
}

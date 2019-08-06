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
package org.comdnmr.fit.gui;

import java.util.ArrayList;
import javafx.collections.ObservableList;
import javafx.scene.Group;
import javafx.scene.Node;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.XYChart;
import javafx.scene.layout.Region;

/**
 *
 * @author Bruce Johnson
 */
public class SSRegion extends Region {

    Group group;
    SSPainter painter;
    NumberAxis xAxis;
    NumberAxis yAxis;
    ArrayList<SecondaryStructure> secondaryStructures;

    public SSRegion() {
        setPrefHeight(30.0);
        setMaxHeight(30.0);
        setMinHeight(30.0);
    }

    public void setChart(XYChart chart, ArrayList<SecondaryStructure> secondaryStructures) {
        xAxis = (NumberAxis) chart.getXAxis();
        yAxis = (NumberAxis) chart.getYAxis();
        this.secondaryStructures = secondaryStructures;
        addSS();
    }

    @Override
    public void layoutChildren() {
        super.layoutChildren();
        addSS();

    }

    public void addSS() {
// 

    }
}

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.cpmgfit2.gui;

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
//        painter = new SSPainter(canvas, xAxis, yAxis, secondaryStructures);
        if (secondaryStructures != null) {
            painter = new SSPainter(xAxis, yAxis, secondaryStructures);
            setPrefHeight(painter.getPrefHeight());
            setMinHeight(painter.getPrefHeight());
            group = null;
            for (Node node : (ObservableList<Node>) getChildren()) {
                if ((node.getId() != null) && node.getId().equals("ssanno")) {
                    group = (Group) node;
                }
            }
            if (group == null) {
                group = new Group();
                group.setId("ssanno");
                getChildren().add(0, group);
            }
            painter.addSS(group);
            //layout();
        }
    }
}

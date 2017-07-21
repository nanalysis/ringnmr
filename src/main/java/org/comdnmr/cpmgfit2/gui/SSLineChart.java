package org.comdnmr.cpmgfit2.gui;

import javafx.scene.chart.NumberAxis;
import java.util.ArrayList;
import javafx.collections.ObservableList;
import javafx.scene.Group;
import javafx.scene.Node;
import javafx.scene.canvas.Canvas;
import javafx.scene.chart.ScatterChart;
import javafx.scene.shape.Polyline;

public class SSLineChart extends ScatterChart {

    Canvas canvas = null;
    SSPainter painter;
    Group group;
    ArrayList<Polyline> polyLines = new ArrayList<>();
    NumberAxis xAxis;
    NumberAxis yAxis;

    public SSLineChart() {
        super(new NumberAxis(), new NumberAxis());
        xAxis = (NumberAxis) getXAxis();
        yAxis = (NumberAxis) getYAxis();
    }

    public void addCanvas() {
        canvas = new Canvas();
        getPlotChildren().add(1, canvas);
    }

    public void addSS(ArrayList<SecondaryStructure> secondaryStructures) {
//        painter = new SSPainter(canvas, xAxis, yAxis, secondaryStructures);
        painter = new SSPainter(xAxis, yAxis, secondaryStructures);
        group = null;
        for (Node node : (ObservableList<Node>) getPlotChildren()) {
            if ((node.getId() != null) && node.getId().equals("ssanno")) {
                group = (Group) node;
            }
        }
        if (group == null) {
            group = new Group();
            group.setId("ssanno");
            getPlotChildren().add(1, group);
        }
        layoutPlotChildren();
    }

    public SSLineChart(NumberAxis xAxis, NumberAxis yAxis) {
        super(xAxis, yAxis);
    }

    @Override
    protected void layoutPlotChildren() {
        //super.getPlotChildren().clear();
        super.layoutPlotChildren();
        if (painter != null) {
            painter.addSS(group);
        }

//        if ((canvas != null) && (painter != null)) {
//            painter.paintSS();
//        }
    }
}

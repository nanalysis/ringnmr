/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.fit.gui;

import javafx.scene.canvas.Canvas;
import javafx.scene.layout.Pane;
import org.nmrfx.chart.XYCanvasChart;

/**
 *
 * @author brucejohnson
 */
public class XYPlotDataPane extends Pane {

    Canvas canvas;
    PlotData chart;

    public XYPlotDataPane() {
        canvas = new Canvas();
        PlotData chart = PlotData.buildChart(canvas);
        getChildren().add(canvas);
        this.chart = chart;
        widthProperty().addListener(e -> updateChart());
        heightProperty().addListener(e -> updateChart());
    }

    public XYCanvasChart getChart() {
        return chart;
    }

    void updateChart() {
        canvas.setWidth(getWidth());
        canvas.setHeight(getHeight());
        chart.drawChart();
    }
}

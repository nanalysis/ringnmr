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
public class XYBarPlotDataPane extends Pane {

    Canvas canvas;

    public XYBarPlotDataPane() {
        canvas = new Canvas();
        getChildren().add(canvas);
        widthProperty().addListener(e -> updateCanvasSize());
        heightProperty().addListener(e -> updateCanvasSize());
    }

    public Canvas getCanvas() {
        return canvas;
    }

    void updateCanvasSize() {
        canvas.setWidth(getWidth());
        canvas.setHeight(getHeight());
    }
}

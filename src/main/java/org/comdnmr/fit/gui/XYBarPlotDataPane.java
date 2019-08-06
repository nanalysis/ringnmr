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

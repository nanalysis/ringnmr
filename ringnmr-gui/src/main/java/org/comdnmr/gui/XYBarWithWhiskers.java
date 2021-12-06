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
package org.comdnmr.gui;

import javafx.scene.Group;
import javafx.scene.control.Tooltip;
import javafx.scene.shape.Line;
import javafx.scene.shape.Rectangle;

/**
 * Candle node used for drawing a candle
 */
public class XYBarWithWhiskers extends Group {

    private Line errorLine = new Line();
    private Rectangle rect = new Rectangle();
    private String seriesStyleClass;
    private String dataStyleClass;
    private Tooltip tooltip = new Tooltip();

    XYBarWithWhiskers(String seriesStyleClass, String dataStyleClass) {
        setAutoSizeChildren(false);
        getChildren().addAll(rect, errorLine);
        this.seriesStyleClass = seriesStyleClass;
        this.dataStyleClass = dataStyleClass;
        updateStyleClasses();
//        tooltip.setGraphic(new TooltipContent());
        Tooltip.install(rect, tooltip);
    }

    public void setSeriesAndDataStyleClasses(String seriesStyleClass, String dataStyleClass) {
        this.seriesStyleClass = seriesStyleClass;
        this.dataStyleClass = dataStyleClass;
        updateStyleClasses();
    }

    public void update(double value, double lowPercentile, double highPercentile, double barWidth) {
        updateStyleClasses();
        errorLine.setStartY(highPercentile);
        errorLine.setEndY(lowPercentile);
        errorLine.setStartX(0 + barWidth / 2);
        errorLine.setEndX(0.0 + barWidth / 2);
//        System.out.println("barwi " + barWidth + " val " + value);
//        rect.setFill(Color.ORANGE);
        rect.setWidth(barWidth);
        rect.setHeight(value);
        updateTooltip(value, lowPercentile, highPercentile);

    }

    public void updateTooltip(double value, double lowPercentile, double highPercentile) {
        tooltip.setText("Value: " + value);
    }

    private void updateStyleClasses() {
        getStyleClass().setAll("xybarchart-bar", seriesStyleClass, dataStyleClass);
        errorLine.getStyleClass().setAll("xybarchart-line", seriesStyleClass, dataStyleClass,
                "close-above-open");
        rect.getStyleClass().setAll("xybarchart-bar", seriesStyleClass, dataStyleClass,
                "open-above-close");
    }
}

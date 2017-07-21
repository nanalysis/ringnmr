package org.comdnmr.cpmgfit2.gui;

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

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.fit.gui;

import java.util.HashSet;
import java.util.List;
import java.util.Optional;
import java.util.Set;
import javafx.collections.ListChangeListener;
import javafx.geometry.Orientation;
import javafx.scene.canvas.Canvas;
import javafx.scene.input.MouseEvent;
import javafx.scene.paint.Color;
import org.comdnmr.data.ResidueInfo;
import org.comdnmr.data.ResidueProperties;
import org.nmrfx.chart.Axis;
import org.nmrfx.chart.XYCanvasBarChart;
import org.nmrfx.chart.XYValue;
import org.nmrfx.graphicsio.GraphicsContextInterface;
import org.nmrfx.graphicsio.GraphicsIOException;

/**
 *
 * @author brucejohnson
 */
public class ResidueChart extends XYCanvasBarChart {

    static Set<String> selectedResidues = new HashSet<>();
    public String currentSeriesName = "";
    ResidueProperties resProps = null;

    public static ResidueChart buildChart(Canvas canvas) {
        Axis xAxis = new Axis(Orientation.HORIZONTAL, 0, 100, 400, 100.0);
        Axis yAxis = new Axis(Orientation.VERTICAL, 0, 100, 100, 400);
        yAxis.setZeroIncluded(true);
        return new ResidueChart(canvas, xAxis, yAxis);
    }

    public ResidueChart(Canvas canvas, final Axis... AXIS) {
        super(canvas, AXIS);
        xAxis = AXIS[0];
        yAxis = AXIS[1];
        init();
    }

    private void init() {
//        canvas.widthProperty().addListener(e -> drawChart());
//        canvas.heightProperty().addListener(e -> drawChart());
        getData().addListener((ListChangeListener) (e -> seriesChanged()));
//        canvas.setOnMouseClicked(e -> mouseClicked(e));
    }

    public boolean mouseClicked(MouseEvent e) {
        Optional<Hit> hitOpt = pickChart(e.getX(), e.getY(), 5);
        boolean hitChart = false;
        if (hitOpt.isPresent()) {
            hitChart = true;
            Hit hit = hitOpt.get();
            System.out.println(hit.toString());
            PyController.mainController.statusBar.setText(hit.toString());
            boolean appendMode = e.isShiftDown();
            XYValue value = hit.getValue();
            int resNum = (int) Math.round(value.getXValue());
            showInfo(hit.getSeries().getName(), hit.getIndex(), resNum, appendMode);
        } else {
            Optional<Integer> intOpt = pickPresenceIndicators(e.getX(), e.getY());
            if (intOpt.isPresent()) {
                hitChart = true;
                boolean appendMode = e.isShiftDown();
                double f = e.getX() / xAxis.getWidth();
                int seriesIndex = (int) Math.floor(getData().size() * f);
                if ((seriesIndex >= 0) && (seriesIndex < getData().size())) {
                    String seriesName = getData().get(seriesIndex).getName();
                    showInfo(seriesName, seriesIndex, intOpt.get(), appendMode);
                }
            }
        }
        return hitChart;
    }

    @Override
    public void annotate(GraphicsContextInterface gC) {
        try {
            drawPresenceIndicators(gC);
        } catch (GraphicsIOException ex) {
        }
    }

    void seriesChanged() {
        selectedResidues.clear();
        drawChart();

    }

    public void setResProps(ResidueProperties resProps) {
        this.resProps = resProps;
    }

    void showInfo(String seriesName, int seriesIndex, int resNum, boolean appendMode) {
        String resString = String.valueOf(resNum);
        if (!appendMode) {
            selectedResidues.clear();
            selectedResidues.add(resString);
        } else if (!selectedResidues.contains(resString)) {
            selectedResidues.add(resString);
        }
        currentSeriesName = seriesName;
        String[] seriesNameParts = seriesName.split("\\|");
        String mapName = seriesNameParts[0];
        ResidueProperties resProps = ChartUtil.residueProperties.get(mapName);
        showInfo(resProps, seriesName);
        drawChart();
    }

    void showInfo() {
        String[] seriesNameParts = currentSeriesName.split("\\|");
        String mapName = seriesNameParts[0];
        ResidueProperties resProps = ChartUtil.residueProperties.get(mapName);
        showInfo(resProps, currentSeriesName);
    }

    void showInfo(ResidueProperties resProps, String seriesName) {
        PyController controller = PyController.mainController;
        PlotData xyCanvasChart = controller.xychart;
        String[] seriesNameParts = seriesName.split("\\|");
        if (seriesNameParts.length < 3) {
            return;
        }
        String mapName = seriesNameParts[0];
        String equationName = seriesNameParts[1];
        String state = seriesNameParts[2];
        state = "*:" + state.substring(2);
        System.out.println("series " + seriesName + " map " + mapName + " eqn " + equationName + " state " + state);
        String[] residues = new String[selectedResidues.size()];
        selectedResidues.toArray(residues);
        controller.showInfo(resProps, equationName, mapName, state, residues, xyCanvasChart);
    }

    Optional<Integer> pickPresenceIndicators(double mouseX, double mouseY) {
        Optional<Integer> result = Optional.empty();
        if (resProps != null) {
            List<ResidueInfo> resValues = resProps.getResidueValues();
            for (ResidueInfo resInfo : resValues) {
                if (resInfo == null) {
                    continue;
                }
                int resNum = resInfo.getResNum();
                double x1 = xAxis.getDisplayPosition(resNum - 0.5) + 1;
                double x2 = xAxis.getDisplayPosition(resNum + 0.5) - 1;
                double y1 = yAxis.getYOrigin() - yAxis.getHeight() + 2;
                double y2 = yAxis.getYOrigin();
                if ((mouseX > x1) && (mouseX < x2)) {
                    if ((mouseY > y1) && (mouseY < y2)) {
                        result = Optional.of(resNum);
                        break;
                    }
                }
            }
        }
        return result;
    }

    void drawPresenceIndicators(GraphicsContextInterface gC) throws GraphicsIOException {
        if (resProps != null) {
            List<ResidueInfo> resValues = resProps.getResidueValues();
            for (ResidueInfo resInfo : resValues) {
                if (resInfo == null) {
                    continue;
                }
                int resNum = resInfo.getResNum();
                double x1 = xAxis.getDisplayPosition(resNum - 0.5) + 1;
                double x2 = xAxis.getDisplayPosition(resNum + 0.5) - 1;
                double y1 = yAxis.getYOrigin() - yAxis.getHeight() + 2;
                double width = x2 - x1;
                double height = yAxis.getHeight() - 4;
                String resName = String.valueOf(resNum);
                if (selectedResidues.contains(resName)) {
                    gC.setFill(Color.ORANGE);
                } else {
                    gC.setFill(Color.LIGHTGRAY);
                }
                gC.fillRect(x1, y1, width, height);
            }
        }

    }
}

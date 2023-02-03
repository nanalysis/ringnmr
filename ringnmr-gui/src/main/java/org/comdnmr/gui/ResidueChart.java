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

import java.util.HashSet;
import java.util.Optional;
import java.util.Set;
import javafx.collections.ListChangeListener;
import javafx.geometry.Orientation;
import javafx.scene.canvas.Canvas;
import javafx.scene.input.MouseEvent;
import javafx.scene.paint.Color;
import org.comdnmr.data.ExperimentSet;
import org.comdnmr.data.RelaxSet;
import org.comdnmr.data.ValueSet;
import org.nmrfx.chart.Axis;
import org.nmrfx.chart.XYCanvasBarChart;
import org.nmrfx.chart.XYValue;
import org.nmrfx.chemistry.MoleculeBase;
import org.nmrfx.chemistry.MoleculeFactory;
import org.nmrfx.chemistry.relax.ResonanceSource;
import org.nmrfx.graphicsio.GraphicsContextInterface;
import org.nmrfx.graphicsio.GraphicsIOException;

/**
 *
 * @author brucejohnson
 */
public class ResidueChart extends XYCanvasBarChart {

    static Set<ResonanceSource> selectedResidues = new HashSet<>();
    public String currentSeriesName = "";
    ValueSet valueSet = null;
    Set<ResonanceSource> dynSources = new HashSet<>();

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
        getData().addListener((ListChangeListener) (e -> seriesChanged()));
    }

    public String getExpSeries(String seriesName) {
        MoleculeBase mol = MoleculeFactory.getActive();

        if (seriesName.startsWith(mol.getName()) && seriesName.contains("_RING")) {
            String testSeries = seriesName.substring(mol.getName().length()+1, seriesName.indexOf("_RING"));
            if (ChartUtil.getResiduePropertyNames().contains(testSeries)) {
                ValueSet valueSet = ChartUtil.getResidueProperty(testSeries);
                seriesName = valueSet.getName() + '|' +
                        "EXPAB" + "|" + "0:0:0" + "|" + "R" + "|" + "0";
            }
        }
        return seriesName;

    }
    public boolean mouseClicked(MouseEvent e) {
        Optional<Hit> hitOpt = pickChart(e.getX(), e.getY(), 5);
        boolean hitChart = false;
        if (hitOpt.isPresent()) {
            hitChart = true;
            Hit hit = hitOpt.get();
            boolean appendMode = e.isShiftDown();
            XYValue value = hit.getValue();
            var series = hit.getSeries();
            Object extraValue = value.getExtraValue();
            ResonanceSource resSource = null;
            if (extraValue instanceof ResonanceSource) {
                resSource = (ResonanceSource) extraValue;
            }
            String seriesName = getExpSeries(series.getName());
            String statusMessage = seriesName + " " + resSource + " " + String.format("%.2f", value.getYValue());
            PyController.mainController.statusBar.setText(statusMessage);
            showInfo(seriesName, hit.getIndex(), resSource, appendMode);
        } else {
            Optional<ResonanceSource> intOpt = pickPresenceIndicators(e.getX(), e.getY());
            if (intOpt.isPresent()) {
                hitChart = true;
                boolean appendMode = e.isShiftDown();
                double f = e.getX() / xAxis.getWidth();
                int seriesIndex = (int) Math.floor(getData().size() * f);
                if ((seriesIndex >= 0) && (seriesIndex < getData().size())) {
                    String seriesName = getData().get(seriesIndex).getName();
                    showInfo(seriesName, seriesIndex, intOpt.get(), appendMode);
                }
            } else {
                selectedResidues.clear();
            }
        }
        return hitChart;
    }

    @Override
    public void annotate(GraphicsContextInterface gC) {
        try {
            if (dynSources != null && !dynSources.isEmpty()) {
                drawPresenceIndicators(gC);
            }
        } catch (GraphicsIOException ex) {
        }
    }

    void seriesChanged() {
        selectedResidues.clear();
        drawChart();

    }

    public void setResProps(ValueSet valueSet) {
        this.valueSet = valueSet;
        dynSources = valueSet.getDynamicsSources();
    }

    public Set<ResonanceSource> getSelectedSources() {
        return selectedResidues;
    }

    void showInfo(String seriesName, int seriesIndex, ResonanceSource resSource, boolean appendMode) {
        if (!appendMode) {
            selectedResidues.clear();
            if (resSource != null) {
                selectedResidues.add(resSource);
            }
        } else if (!selectedResidues.contains(resSource)) {
            if (resSource != null) {
                selectedResidues.add(resSource);
            }
        }
        System.out.println(" Res source is " + resSource);
        currentSeriesName = seriesName;
        String[] seriesNameParts = seriesName.split("\\|");
        String mapName = seriesNameParts[0];
        valueSet = ChartUtil.getResidueProperty(mapName);
        dynSources = valueSet.getDynamicsSources();
        showInfo(seriesName);
        drawChart();
    }

    void showInfo() {
        String[] seriesNameParts = currentSeriesName.split("\\|");
        String mapName = seriesNameParts[0];
        valueSet = ChartUtil.getResidueProperty(mapName);
        dynSources = valueSet.getDynamicsSources();
        showInfo(currentSeriesName);
    }

    void showInfo(String seriesName) {
        PyController controller = PyController.mainController;
        PlotData xyCanvasChart = controller.xychart;
        String[] seriesNameParts = seriesName.split("\\|");
        System.out.println("show info se " + seriesName + " " + seriesNameParts.length + " valueset " + valueSet);
        if (seriesNameParts.length < 3) {
            return;
        }
        String mapName = seriesNameParts[0];
        String equationName = seriesNameParts[1];
        String state = seriesNameParts[2];
        state = "*:" + state.substring(2);
//        System.out.println("series " + seriesName + " map " + mapName + " eqn " + equationName + " state " + state);
        ResonanceSource[] resSources = new ResonanceSource[selectedResidues.size()];
        selectedResidues.toArray(resSources);
        controller.chartInfo.valueSet = valueSet;
        controller.chartInfo.setResidues(resSources);
        controller.chartInfo.state = state;
        controller.chartInfo.mapName = mapName;
        controller.chartInfo.equationName = equationName;

        if (valueSet instanceof ExperimentSet) {
            if (!selectedResidues.isEmpty()) {
                ExperimentSet experimentSet = (ExperimentSet) valueSet;
                controller.chartInfo.experimentalResult = experimentSet.getExperimentResult(resSources[0]);
            }
        } else if (valueSet instanceof RelaxSet) {
        }
        controller.showInfo(controller.chartInfo, xyCanvasChart);
    }

    Optional<ResonanceSource> pickPresenceIndicators(double mouseX, double mouseY) {
        Optional<ResonanceSource> result = Optional.empty();
        if (dynSources != null) {
            for (var dynSource : dynSources) {
                int resNum = dynSource.getAtom().getResidueNumber();
                double x1 = xAxis.getDisplayPosition(resNum - 0.5) + 1;
                double x2 = xAxis.getDisplayPosition(resNum + 0.5) - 1;
                double y1 = yAxis.getYOrigin() - yAxis.getHeight() + 2;
                double y2 = yAxis.getYOrigin();
                if ((mouseX > x1) && (mouseX < x2)) {
                    if ((mouseY > y1) && (mouseY < y2)) {
                        result = Optional.of(dynSource);
                        break;
                    }
                }
            }
        }
        return result;
    }

    void drawPresenceIndicators(GraphicsContextInterface gC) throws GraphicsIOException {
        if (dynSources != null) {
            for (var dynSource : dynSources) {
                int resNum = dynSource.getAtom().getResidueNumber();
                double x1 = xAxis.getDisplayPosition(resNum - 0.5) + 1;
                double x2 = xAxis.getDisplayPosition(resNum + 0.5) - 1;
                double y1 = yAxis.getYOrigin() - yAxis.getHeight() + 2;
                double width = x2 - x1;
                double height = yAxis.getHeight() - 4;
                if (selectedResidues.contains(dynSource)) {
                    gC.setFill(Color.ORANGE);
                } else {
                    gC.setFill(Color.LIGHTGRAY);
                }
                gC.fillRect(x1, y1, width, height);
            }
        }
    }
}

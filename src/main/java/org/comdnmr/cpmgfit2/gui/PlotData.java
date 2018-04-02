package org.comdnmr.cpmgfit2.gui;

import java.io.BufferedReader;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.scene.canvas.Canvas;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.XYChart;
import javafx.scene.chart.XYChart.Series;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Set;
import javafx.geometry.Side;
import javafx.scene.Node;
import javafx.scene.chart.ScatterChart;
import javafx.scene.chart.XYChart.Data;
import javafx.scene.shape.Polyline;
import javafx.scene.shape.Circle;
import javafx.scene.shape.Line;
import javafx.scene.Group;
import javafx.scene.control.Label;

import javafx.scene.paint.Color;
import org.comdnmr.cpmgfit2.calc.PlotEquation;
import org.comdnmr.cpmgfit2.calc.ResidueData;

public class PlotData extends ScatterChart {

    NumberAxis xAxis;
    NumberAxis yAxis;
    ObservableList<XYChart.Series<Double, Double>> simData = FXCollections.observableArrayList();
    String fileName;
    ArrayList<Polyline> polyLines = new ArrayList<>();
    ArrayList<PlotEquation> plotEquations = new ArrayList<>();
    static Color[] colors = {Color.ORANGE, Color.BLUE, Color.RED, Color.GREEN, Color.GRAY};

    public PlotData(NumberAxis xAxis, NumberAxis yAxis) {
        super(xAxis, yAxis);
        this.xAxis = xAxis;
        this.yAxis = yAxis;
        init();
    }

    public PlotData() {
        super(new NumberAxis(), new NumberAxis());
        xAxis = (NumberAxis) getXAxis();
        yAxis = (NumberAxis) getYAxis();
        init();
    }

    void init() {
        xAxis.setAutoRanging(true);
        xAxis.setUpperBound(1000.0);
        yAxis.setUpperBound(60.0);
        xAxis.setLabel("\u03BD(cpmg)");
        yAxis.setLabel("R2(\u03BD)");
        xAxis.setLabel("Time (s)");
        yAxis.setLabel("Intensity");
        setTitle("CurveFit");
        setPrefHeight(200);
        setLegendSide(Side.BOTTOM);
        setLegendVisible(true);
        xAxis.setAnimated(false);
        yAxis.setAnimated(false);
        setAnimated(false);
        setNodeListeners(this);
    }

    public void addCanvas(Canvas canvas) {
        getPlotChildren().add(1, canvas);
    }

    public void setEquations(List<PlotEquation> plotEquations) {
        this.plotEquations.clear();
        this.plotEquations.addAll(plotEquations);
    }

    @Override
    protected void layoutPlotChildren() {
        if (true || (polyLines.size() != plotEquations.size())) {
            getPlotChildren().removeAll(polyLines);
            polyLines.clear();
            int nLines = plotEquations.size();
            for (int i = 0; i < nLines; i++) {
                Polyline polyLine = new Polyline();
                polyLine.setStroke(colors[Math.min(colors.length - 1, i)]);
                polyLines.add(polyLine);
                getPlotChildren().add(0, polyLine);
            }

        }
        super.layoutPlotChildren();
        paintLines(polyLines, ((NumberAxis) getXAxis()), ((NumberAxis) getYAxis()));
        setNodeListeners(this);
    }

    void setNodeListeners(XYChart chart) {
        ObservableList<XYChart.Series<Double, Double>> data = chart.getData();
        for (XYChart.Series<Double, Double> series : data) {
            int j = 0;
            for (Data<Double, Double> item : series.getData()) {
                Node node = item.getNode();
                if (node != null) {
//                    node.onMouseClickedProperty().addListener(e -> dumpNode(item));
                    node.setOnMouseClicked(e -> dumpNode(series.getName(), item));
                    if (node instanceof Group) {
                        Group group = (Group) node;
                        Line line = (Line) group.getChildren().get(1);
                        double errorY = getScaledError(item);
                        line.setStartY(-errorY);
                        line.setEndY(errorY);
                    }
                }
            }

        }
    }

    void dumpNode(String seriesName, Data<Double, Double> data) {
        Object extraValue = data.getExtraValue();
        if (extraValue instanceof ResidueData.DataValue) {
            ResidueData.DataValue dataValue = (ResidueData.DataValue) extraValue;
            PyController.mainController.selectTableRow(seriesName, dataValue.getIndex());

        }

    }

    void paintLinesSeries(ArrayList<Polyline> polyLines, NumberAxis xAxis, NumberAxis yAxis) {
        int i = 0;
        for (XYChart.Series<Double, Double> series : simData) {
            ArrayList<Double> points = new ArrayList<>();
            Polyline polyLine = polyLines.get(i);
            polyLine.getPoints().clear();
            for (Data<Double, Double> xyData : series.getData()) {
                double x = xAxis.getDisplayPosition(xyData.getXValue());
                double y = yAxis.getDisplayPosition(xyData.getYValue());
                points.add(x);
                points.add(y);
            }
            polyLine.getPoints().addAll(points);
            i++;
        }
    }

    void paintLines(ArrayList<Polyline> polyLines, NumberAxis xAxis, NumberAxis yAxis) {
        //double[] fields = {11.7,14.04};
        //double[] par = {6.315686252794991, 8.612586767595255, 0.7538704344970752};
        //equation.setFieldRef(fields[0]);
        int iLine = 0;
        if (polyLines.isEmpty()) {
            return;
        }
        double fieldRef = 1.0;
        for (PlotEquation plotEquation : plotEquations) {
            if (plotEquation == null) {
                continue;
            }
            ArrayList<Double> points = new ArrayList<>();
            Polyline polyLine = polyLines.get(iLine);
            polyLine.getPoints().clear();
            double min = xAxis.getLowerBound();
            double max = xAxis.getUpperBound();
            int nIncr = 100;
            double delta = (max - min) / nIncr;
            if (iLine == 0) {
                fieldRef = plotEquation.getExtra(0);
            }
//            plotEquation.equation.setFieldRef(fieldRef);
            for (int i = 1; i < nIncr - 1; i++) {
                double xValue = min + i * delta;
                double x = xAxis.getDisplayPosition(xValue);
                double yValue = plotEquation.calculate(xValue, plotEquation.getExtra(0) / fieldRef);
                double y = yAxis.getDisplayPosition(yValue);
                points.add(x);
                points.add(y);
            }
            polyLine.getPoints().addAll(points);
            iLine++;

        }

    }

    private ObservableList<XYChart.Series<Double, Double>> loadChartData(String[] residues) throws IOException {
        ObservableList<XYChart.Series<Double, Double>> data = FXCollections.observableArrayList();
        Path path = Paths.get(".", fileName);
        Series<Double, Double> series = null;
        simData.clear();
        boolean inData = false;
        double[] xValues = null;
        double[] yValues = null;
        try (BufferedReader fileReader = Files.newBufferedReader(path)) {
            String line = fileReader.readLine();
            double field = Double.parseDouble(line.trim());
            line = fileReader.readLine();
            double time = Double.parseDouble(line.trim());
            while (true) {
                line = fileReader.readLine();
                if (line == null) {
                    break;
                }
                String sline = line.trim();
                String[] sfields = sline.split("\t");
                if (line.startsWith("#")) {
                } else if (line.startsWith("vCPMG")) {
                    xValues = new double[sfields.length - 2];
                    for (int i = 2; i < sfields.length; i++) {
                        xValues[i - 2] = Double.parseDouble(sfields[i]);
                    }
                } else if (Arrays.stream(residues).anyMatch(e -> e.equals(sfields[0]))) {
                    series = new Series<>();
                    data.add(series);

                    series.setName(sfields[0]);
                    double ref = Double.parseDouble(sfields[1]);
                    yValues = new double[sfields.length - 2];
                    for (int i = 2; i < sfields.length; i++) {
                        yValues[i - 2] = -Math.log(Double.parseDouble(sfields[i]) / ref) / time;
                        //-math.log(v/ref)/tCPMG
                    }
                    for (int i = 0; i < xValues.length; i++) {
                        if (series != null) {
                            double x = xValues[i];
                            double y = yValues[i];
                            XYChart.Data dataPoint = new XYChart.Data(x, y);
                            series.getData().add(dataPoint);
                        }

                    }
                }
            }
        }
        return data;
    }

    void doLine(Series<Integer, Double> series, String line, double scale) {
        String[] fields = line.split(" ");
        Integer x = Integer.parseInt(fields[0]);
        Double y = Double.parseDouble(fields[1]) * scale;
        series.getData().add(new XYChart.Data(x, y));
    }

    private double getScaledError(XYChart.Data item) {
        Object extraValue = item.getExtraValue();
        double errorY = 0.0;
        if (extraValue instanceof ResidueData.DataValue) {
            ResidueData.DataValue dataValue = (ResidueData.DataValue) extraValue;
            double error = dataValue.getError();
            double yValue = dataValue.getY();
            double errorY2 = yAxis.getDisplayPosition(yValue + error);
            double errorY1 = yAxis.getDisplayPosition(yValue - error);
            errorY = Math.abs(errorY2 - errorY1) / 2.0;
        }
        return errorY;

    }

    private Node createNode(XYChart.Data item, int seriesIndex) {
        Object extraValue = item.getExtraValue();
        double errorY = getScaledError(item);
        Group g = new Group();
        Circle circle = new Circle(4.0);
        circle.setFill(colors[Math.min(colors.length - 1, seriesIndex)]);
        Line line = new Line(0, -errorY, 0, errorY);
        g.getChildren().add(circle);
        g.getChildren().add(line);
        return g;
    }

    @Override
    protected void seriesAdded(XYChart.Series series, int seriesIndex) {
        super.seriesAdded(series, seriesIndex);
        for (int j = 0; j < series.getData().size(); j++) {
            XYChart.Data item = (XYChart.Data) series.getData().get(j);
            Node node = createNode(item, seriesIndex);
            item.setNode(node);
            super.dataItemAdded(series, j, item);
        }
        updateLegend();
    }

    protected void updateLegend() {
        super.updateLegend();
        Set<Node> items = lookupAll(".chart-legend-item");
        int it = 0;
        for (Node item : items) {
            Label label = (Label) item;
            Circle circle = new Circle(4.0);
            circle.setFill(colors[Math.min(colors.length - 1, it)]);
            label.setGraphic(circle);
            it++;
        }
    }

//    @Override
//    protected void seriesAdded(XYChart.Series series, int seriesIndex) {
//        // handle any data already in series
//        for (int j = 0; j < series.getData().size(); j++) {
//            XYChart.Data item = (XYChart.Data) series.getData().get(j);
////            Circle circle = new Circle(20.0);
////            circle.setFill(Color.ORANGE);
////            item.setNode(circle);
//            item.getNode().setOnMouseEntered(new EventHandler<MouseEvent>() {
//                @Override
//                public void handle(MouseEvent mouseEvent) {
//                    System.out.println("enter");
//                }
//            });
//            item.getNode().setOnMouseExited(new EventHandler<MouseEvent>() {
//                @Override
//                public void handle(MouseEvent mouseEvent) {
//                    System.out.println("exit");
//                }
//            });
//
//        }
//    }
}

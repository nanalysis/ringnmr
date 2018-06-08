package org.comdnmr.fit.gui;

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
import java.util.HashMap;
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
import org.comdnmr.fit.calc.ExperimentData;
import org.comdnmr.fit.calc.PlotEquation;
import org.comdnmr.fit.calc.ResidueData;
import org.comdnmr.fit.calc.ResidueProperties;
import static org.comdnmr.fit.gui.ChartUtil.residueProperties;

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
        //xAxis.setAutoRanging(true);
        xAxis.setUpperBound(1000.0);
        yAxis.setUpperBound(60.0);
        xAxis.setLabel("Time (s)");
        yAxis.setLabel("Intensity");
        xAxis.setLabel("\u03BD (cpmg)");
        yAxis.setLabel("R2 (\u03BD)");
        setTitle("CPMG");
        setPrefHeight(200);
        setLegendSide(Side.BOTTOM);
        setLegendVisible(true);
        xAxis.setAnimated(false);
        yAxis.setAnimated(false);
        setAnimated(false);
        setNodeListeners(this);
        setHorizontalZeroLineVisible(false);
        setVerticalZeroLineVisible(false);
    }

    public void setNames(String title, String xName, String yName, String yPad) {
        setTitle(title);
        xAxis.setLabel(xName);
        yAxis.setLabel(yName);
        getYAxis().lookup(".axis-label").setStyle("-fx-label-padding: 0 0 " + yPad + " 0;");
    }

    public void setBounds(double xLower, double xUpper, double yLower, double yUpper, double xtick, double ytick) {
        xAxis.setAutoRanging(false);
        yAxis.setAutoRanging(false);
        xAxis.setLowerBound(xLower);
        xAxis.setUpperBound(xUpper);
        yAxis.setLowerBound(yLower);
        yAxis.setUpperBound(yUpper);
        xAxis.setTickUnit(xtick);
        yAxis.setTickUnit(ytick);
    }

    public void addCanvas(Canvas canvas) {
        getPlotChildren().add(1, canvas);
    }

    public int getNumPlots() {
        int numPlots = getData().size();
        return numPlots;
    }

    public int getNumEquations() {
        int numEquations = plotEquations.size();
        return numEquations;
    }

    public void setEquations(List<PlotEquation> plotEquations) {
        this.plotEquations.clear();
        this.plotEquations.addAll(plotEquations);
    }

    public void clear() {
        plotEquations.clear();
        getData().clear();
    }

    @Override
    protected void layoutPlotChildren() {
        if (true || (polyLines.size() != plotEquations.size())) {
            getPlotChildren().removeAll(polyLines);
            polyLines.clear();
            int nLines = plotEquations.size();
            //System.out.println("layoutPlotChildren No. Eqns = " + nLines);
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
            //System.out.println("paintLines No. Eqns = " + plotEquations.size());
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
            double[] extras = plotEquation.getExtras();
            double[] ax = new double[extras.length];
//System.out.println("extras " + extras.length);
            for (int i = 1; i < nIncr - 1; i++) {
                double xValue = min + i * delta;
                double x = xAxis.getDisplayPosition(xValue);
                ax[0] = xValue;
                for (int j = 1; j < extras.length; j++) {
                    ax[j] = extras[j];
                }
                double yValue = plotEquation.calculate(ax, plotEquation.getExtra(0) / fieldRef);
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
    // Returns a list of maps that store data about each plot on cpmgfit

    public HashMap<String, Object>[] getPlottedData() {
        int numPlots = getNumPlots();
        int numFits = getNumEquations();

        if (numPlots == numFits) {
            // The fitted curve belongs to the raw data collected
            HashMap<String, Object>[] plottedData = new HashMap[numPlots];
            for (int i = 0; i < numPlots; i++) {
                HashMap<String, Object> plotDatum = new HashMap<>();
                plotDatum.put("rawData", returnLine(i));
                plotDatum.put("graphTitle", getSeriesName(i));
                plotDatum.put("fittedData", returnEquation(i));
                plottedData[i] = plotDatum;
            }
            return plottedData;

        }
        //Will break if numPlots != numFits
        return null;
    }

    // Returns a map of parameters about the graph
    public HashMap<String, Object> getGraphData() {
        HashMap<String, Object> graphData = new HashMap<>();

        String title = getTitle();
        String xLabel = xAxis.getLabel();
        String yLabel = yAxis.getLabel();
        double xMin = xAxis.lowerBoundProperty().getValue();
        double xMax = xAxis.upperBoundProperty().getValue();
        double yMin = yAxis.lowerBoundProperty().getValue();
        double yMax = yAxis.upperBoundProperty().getValue();

        List<Double> ranges = new ArrayList(Arrays.asList(xMin, xMax, yMin, yMax));
        graphData.put("title", title);
        graphData.put("xlabel", xLabel);
        graphData.put("ylabel", yLabel);
        graphData.put("ranges", ranges);
        graphData.put("colors", PlotData.colors);
        return graphData;
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
        ResidueProperties residueProps = residueProperties.get("cest");
        ExperimentData expData = null;
        if (residueProps != null) {
            expData = residueProps.getExperimentData("cest"); // fixme
        }
        for (int j = 0; j < series.getData().size(); j++) {
            XYChart.Data item = (XYChart.Data) series.getData().get(j);
            ResidueData.DataValue dataValue = (ResidueData.DataValue) item.getExtraValue();
            double x1 = dataValue.getX1();
            Node node = createNode(item, seriesIndex);
            double x1val = Math.round(100.0 * x1 / (2 * Math.PI)) / 100.0;
            if (expData != null && expData.getExtras().contains(x1val)) {
                node = createNode(item, expData.getExtras().indexOf(x1val));
            }
            item.setNode(node);
            super.dataItemAdded(series, j, item);
        }
        updateLegend();
//        if (expData != null) {
//                cestLegend(expData.getExtras());
//            }
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

//    protected void cestLegend(List<Double> b1fields) {
//        super.updateLegend();
//        Set<Node> items = b1fields;
//        int it = 0;
//        for (Node item : items) {
//            Label label = (Label) item;
//            Circle circle = new Circle(4.0);
//            circle.setFill(colors[Math.min(colors.length - 1, it)]);
//            label.setGraphic(circle);
//            it++;
//        }
//    }
    public String getSeriesName(int iSeries) {
        ObservableList<XYChart.Series<Double, Double>> data = getData();
        XYChart.Series<Double, Double> series = data.get(iSeries);
        String seriesName = series.getName();  // Can get residue by splitting semicolon
        return seriesName;
    }

    public void dumpLine(int iSeries) {
        ArrayList<ArrayList<Double>> lineData = returnLine(iSeries);
        lineData.forEach((pointData -> {
            String outputLine = String.format("%14.4g %14.4g %14.4g", pointData.get(0), pointData.get(1), pointData.get(2));
            System.out.println(pointData);
        }));
    }

    public ArrayList<ArrayList<Double>> returnLine(int iSeries) {
        ObservableList<XYChart.Series<Double, Double>> data = getData();
        XYChart.Series<Double, Double> series = data.get(iSeries);
        ArrayList<ArrayList<Double>> lineData = new ArrayList<>(series.getData().size());
        series.getData().forEach((value) -> {
            Double x = value.getXValue();
            Double y = value.getYValue();
            Object extraValue = value.getExtraValue();
            if (extraValue instanceof ResidueData.DataValue) {
                Double error = ((ResidueData.DataValue) extraValue).getError();
                ArrayList<Double> pointData = new ArrayList<>(3);
                pointData.add(x);
                pointData.add(y);
                pointData.add(error);
                lineData.add(pointData);
            }
        });
        return lineData;
    }

    public void dumpEquation(int iLine) {
        ArrayList<ArrayList<Double>> equationData = returnEquation(iLine);
        equationData.forEach((pointData -> {
            String outputLine = String.format("%14.4g %14.4g %14.4g", pointData.get(0), pointData.get(1), pointData.get(2));
            System.out.println(pointData);
        }));
    }

    public ArrayList<ArrayList<Double>> returnEquation(int iLine) {
        PlotEquation plotEquation = plotEquations.get(iLine);
        double fieldRef;
        if (iLine == 0) {
            fieldRef = plotEquation.getExtra(0);
        } else { //Prevents another call if iLine = 0; probably an easier way to do this
            PlotEquation refEquation = plotEquations.get(0);
            fieldRef = refEquation.getExtra(0);
        }
        double min = xAxis.getLowerBound();
        double max = xAxis.getUpperBound();
        int nIncr = 100;
        double delta = (max - min) / nIncr;
        ArrayList<ArrayList<Double>> equationData = new ArrayList<>(nIncr - 1);
        double[] ax = new double[1];
        for (int i = 1; i < nIncr - 1; i++) {
            ArrayList<Double> pointData = new ArrayList<>(3);
            double xValue = min + i * delta;
            ax[0] = xValue;
            double yValue = plotEquation.calculate(ax, plotEquation.getExtra(0) / fieldRef);
            pointData.add(xValue);
            pointData.add(yValue);
            equationData.add(pointData);
        }
        return equationData;
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

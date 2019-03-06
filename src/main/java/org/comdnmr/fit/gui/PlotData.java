package org.comdnmr.fit.gui;

import java.io.BufferedReader;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.scene.canvas.Canvas;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Optional;
import javafx.collections.ListChangeListener;
import javafx.geometry.Orientation;
import javafx.scene.canvas.GraphicsContext;
import javafx.scene.input.MouseEvent;

import javafx.scene.paint.Color;
import org.comdnmr.fit.calc.ExperimentData;
import org.comdnmr.fit.calc.PlotEquation;
import org.comdnmr.fit.calc.ResidueData;
import org.comdnmr.fit.calc.ResidueProperties;
import static org.comdnmr.fit.gui.ChartUtil.residueProperties;
import org.nmrfx.chart.Axis;
import org.nmrfx.chart.DataSeries;

import org.nmrfx.chart.XYCanvasChart;
import org.nmrfx.chart.XYValue;
import org.nmrfx.graphicsio.GraphicsContextInterface;
import org.nmrfx.graphicsio.GraphicsContextProxy;
import org.nmrfx.graphicsio.GraphicsIOException;
import org.nmrfx.graphicsio.SVGGraphicsContext;

public class PlotData extends XYCanvasChart {

    ObservableList<DataSeries> simData = FXCollections.observableArrayList();
    String fileName;
    ObservableList<GUIPlotEquation> plotEquations = FXCollections.observableArrayList();

    public static final Color[] colors = {
        Color.web("#1b9e77"),
        Color.web("#d95f02"),
        Color.web("#7570b3"),
        Color.web("#e7298a"),
        Color.web("#66a61e"),
        Color.web("#e6ab02"),
        Color.web("#a6761d"),
        Color.web("#666666"),
        Color.web("#ff7f00"),
        Color.web("#6a3d9a"),};

    public PlotData(Canvas canvas, final Axis... AXIS) {
        super(canvas, AXIS);
        init();
    }

    public static PlotData buildChart(Canvas canvas) {
        Axis xAxis = new Axis(Orientation.HORIZONTAL, 0, 100, 400, 100.0);
        Axis yAxis = new Axis(Orientation.VERTICAL, 0, 100, 100, 400);
        return new PlotData(canvas, xAxis, yAxis);
    }

    void init() {
        //xAxis.setAutoRanging(true);
        //yAxis.setAutoRanging(true);
        xAxis.setUpperBound(1000.0);
        yAxis.setUpperBound(60.0);
        yAxis.setZeroIncluded(true);
        xAxis.setLabel("Time (s)");
        yAxis.setLabel("Intensity");
        xAxis.setLabel("\u03BD (cpmg)");
        yAxis.setLabel("R2 (\u03BD)");
        plotEquations.addListener((ListChangeListener) (e -> drawChart()));
        getCanvas().setOnMouseClicked(e -> mouseClicked(e));
//        setTitle("CPMG");
//        setPrefHeight(200);
//        setLegendSide(Side.BOTTOM);
//        setLegendVisible(true);
//        xAxis.setAnimated(false);
//        yAxis.setAnimated(false);
//        setAnimated(false);
//        setNodeListeners(this);
//        setHorizontalZeroLineVisible(false);
//        setVerticalZeroLineVisible(false);
    }

    public void setNames(String title, String xName, String yName, String yPad) {
        setTitle(title);
        xAxis.setLabel(xName);
        yAxis.setLabel(yName);
    }

    public void setBounds(double xLower, double xUpper, double yLower, double yUpper, double xtick, double ytick) {
//        xAxis.setAutoRanging(false);
//        yAxis.setAutoRanging(false);
        xAxis.setLowerBound(xLower);
        xAxis.setUpperBound(xUpper);
        yAxis.setLowerBound(yLower);
        yAxis.setUpperBound(yUpper);
        drawChart();
//        xAxis.setTickUnit(xtick);
//        yAxis.setTickUnit(ytick);
    }

    void mouseClicked(MouseEvent e) {
        Optional<Hit> hitOpt = pickChart(e.getX(), e.getY(), 5);
        if (hitOpt.isPresent()) {
            Hit hit = hitOpt.get();
            PyController.mainController.selectTableRow(hit.getSeries().getName(), hit.getIndex());
            PyController.mainController.statusBar.setText(hit.toString());

        }
    }

    public int getNumPlots() {
        int numPlots = getData().size();
        return numPlots;
    }

    public int getNumEquations() {
        int numEquations = plotEquations.size();
        return numEquations;
    }

    public void setEquations(List<GUIPlotEquation> plotEquations) {
        this.plotEquations.clear();
        this.plotEquations.addAll(plotEquations);
        drawChart();
    }

    public void clear() {
        plotEquations.clear();
        getData().clear();
    }

    @Override
    public void drawChart() {
        super.drawChart();
        try {
            paintLines();
        } catch (GraphicsIOException ex) {
            ex.printStackTrace();
        }
    }

    protected void exportVectorGraphics(SVGGraphicsContext svgGC) throws GraphicsIOException {
        svgGC.save();
        svgGC.clip();
        svgGC.beginPath();
        drawChart(svgGC);
        paintLines(svgGC);
        svgGC.restore();
    }

    void dumpNode(String seriesName, XYValue value) {
        Object extraValue = value.getExtraValue();
        if (extraValue instanceof ResidueData.DataValue) {
            ResidueData.DataValue dataValue = (ResidueData.DataValue) extraValue;
            PyController.mainController.selectTableRow(seriesName, dataValue.getIndex());
        }
    }

    public double[] getXBounds() {
        double[] bounds = {xAxis.getLowerBound(), xAxis.getUpperBound()};
        return bounds;
    }

    public double[] getYBounds() {
        double[] bounds = {yAxis.getLowerBound(), yAxis.getUpperBound()};
        return bounds;
    }

    void paintLines() throws GraphicsIOException {
        GraphicsContext gCC = getCanvas().getGraphicsContext2D();
        GraphicsContextInterface gC = new GraphicsContextProxy(gCC);
        paintLines(gC);
    }

    void paintLines(GraphicsContextInterface gC) throws GraphicsIOException {

        int nIncr = 256;
        double[] xValues = new double[nIncr];
        double[] yValues = new double[nIncr];
        for (GUIPlotEquation plotEquation : plotEquations) {
            if (plotEquation == null) {
                continue;
            }
            double min = Math.max(xAxis.getLowerBound(), plotEquation.getMinX());
            double max = Math.min(xAxis.getUpperBound(), plotEquation.getMaxX());

            double delta = (max - min) / (nIncr + 1);
            double[] extras = plotEquation.getExtras();
            double[] ax = new double[extras.length];
            for (int i = 0; i < nIncr; i++) {
                double xValue = min + (i + 1) * delta;
                double x = xAxis.getDisplayPosition(xValue);
                ax[0] = xValue;
                for (int j = 1; j < extras.length; j++) {
                    ax[j] = extras[j];
                }
                double yValue = plotEquation.calculate(ax, plotEquation.getExtra(0));
               // yValue /= plotEquation.getScaleValue() / 100.0;
                double y = yAxis.getDisplayPosition(yValue);
                xValues[i] = x;
                yValues[i] = y;
            }
            gC.setStroke(plotEquation.getColor());
            gC.strokePolyline(xValues, yValues, nIncr);
        }
    }

    private ObservableList<DataSeries> loadChartData(String[] residues) throws IOException {
        ObservableList<DataSeries> data = FXCollections.observableArrayList();
        Path path = Paths.get(".", fileName);
        DataSeries series = null;
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
                    series = new DataSeries();
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
                            XYValue dataPoint = new XYValue(x, y);
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

//    void doLine(Series<Integer, Double> series, String line, double scale) {
//        String[] fields = line.split(" ");
//        Integer x = Integer.parseInt(fields[0]);
//        Double y = Double.parseDouble(fields[1]) * scale;
//        series.getData().add(new XYChart.Data(x, y));
//    }
    private double getScaledError(XYValue item) {
        Object extraValue = item.getExtraValue();
        double errorY = 0.0;
        if (extraValue instanceof ResidueData.DataValue) {
            ResidueData.DataValue dataValue = (ResidueData.DataValue) extraValue;
            double error = dataValue.getError();
            double yValue = dataValue.getY();
            double errorY2 = yAxis.getDisplayPosition(yValue + error);
            double errorY1 = yAxis.getDisplayPosition(yValue - error);
            errorY = Math.abs(errorY2 - errorY1) / 2.0;
        } else if (extraValue instanceof Double) {
            double yValue = ((Double) item.getYValue());
            double error = (Double) extraValue;
            double errorY2 = yAxis.getDisplayPosition(yValue + error);
            double errorY1 = yAxis.getDisplayPosition(yValue - error);
            errorY = Math.abs(errorY2 - errorY1) / 2.0;
        }
        return errorY;

    }

    protected void seriesAdded(DataSeries series, int seriesIndex) {
        ResidueProperties residueProps = residueProperties.get("cest");
        ExperimentData expData = null;
        if (residueProps != null) {
            expData = residueProps.getExperimentData("cest"); // fixme
        }
        for (int j = 0; j < series.getData().size(); j++) {
            XYValue item = (XYValue) series.getData().get(j);
        }
//        if (expData != null) {
//                cestLegend(expData.getExtras());
//            }
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
        ObservableList<DataSeries> data = getData();
        DataSeries series = data.get(iSeries);
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
        ObservableList<DataSeries> data = getData();
        DataSeries series = data.get(iSeries);
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

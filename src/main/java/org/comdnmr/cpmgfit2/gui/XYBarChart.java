package org.comdnmr.cpmgfit2.gui;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import javafx.animation.FadeTransition;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.event.ActionEvent;
import javafx.scene.Group;
import javafx.scene.Node;
import javafx.scene.chart.Axis;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.XYChart;
import javafx.scene.input.ContextMenuEvent;
import javafx.scene.input.MouseEvent;
import javafx.scene.paint.Color;
import javafx.scene.shape.Path;
import javafx.scene.shape.Rectangle;
import javafx.util.Duration;
import org.comdnmr.cpmgfit2.calc.ExperimentData;
import org.comdnmr.cpmgfit2.calc.PlotEquation;
import org.comdnmr.cpmgfit2.calc.ResidueData;
import org.comdnmr.cpmgfit2.calc.ResidueInfo;
import org.comdnmr.cpmgfit2.calc.ResidueProperties;

public class XYBarChart extends XYChart<Number, Number> {

    NumberAxis xAxis;
    NumberAxis yAxis;
    SSPainter painter;
    Group group = null;
    Group presenceGroup = null;
    boolean unifyYAxes = true;
    List<String> selectedResidues = new ArrayList<>();
    String currentSeriesName = "";
    ResidueProperties resProps = null;

    // -------------- CONSTRUCTORS ----------------------------------------------
    public XYBarChart() {
        super(new NumberAxis(), new NumberAxis());
        xAxis = (NumberAxis) getXAxis();
        yAxis = (NumberAxis) getYAxis();
//        xAxis.setAutoRanging(false);
        this.xAxis.setForceZeroInRange(false);

        Object rescss = XYBarChart.class.getResource("XYBarChart.css");
        getStylesheets().add(XYBarChart.class.getResource("XYBarChart.css").toExternalForm());

//        XYChart.Series<Number, Number> series = new XYChart.Series<Number, Number>();
        ObservableList<XYChart.Series<Number, Number>> data = FXCollections.observableArrayList();
        setData(data);
        setPrefHeight(200);
        this.setOnMousePressed(e -> mousePressed(e));
        yAxis.setOnContextMenuRequested(e -> showAxisMenu(e));

//        Insets insets = new Insets(0, 0, 100, 0);
//        setPadding(insets);
    }

    void showAxisMenu(ContextMenuEvent event) {
        Node mouseNode = (Node) event.getSource();
        // PyController.mainController.axisMenu.show(mouseNode.getScene().getWindow(), event.getScreenX(), event.getScreenY());
    }

    /**
     * Construct a new CandleStickChart with the given axis.
     *
     * @param xAxis The x axis to use
     * @param yAxis The y axis to use
     */
    public XYBarChart(Axis<Number> xAxis, Axis<Number> yAxis) {
        super(xAxis, yAxis);
        Object rescss = XYBarChart.class.getResource("XYBarChart.css");
        getStylesheets().add(XYBarChart.class.getResource("XYBarChart.css").toExternalForm());
        setAnimated(false);
        xAxis.setAnimated(false);
        yAxis.setAnimated(false);
        double xMin = 0;
        double xMax = 100.0;
        this.xAxis.setLowerBound(xMin);
        this.xAxis.setUpperBound(xMax);
        this.xAxis.setTickUnit(10);
        this.xAxis.setMinorTickCount(5);
        this.yAxis.setMinWidth(75.0);
        this.yAxis.setPrefWidth(75.0);

//        xAxis.setAutoRanging(true);
//        this.xAxis.setForceZeroInRange(false);
    }

    /**
     * Construct a new CandleStickChart with the given axis and data.
     *
     * @param xAxis The x axis to use
     * @param yAxis The y axis to use
     * @param data The data to use, this is the actual list used so any changes to it will be reflected in the chart
     */
    public XYBarChart(Axis<Number> xAxis, Axis<Number> yAxis, ObservableList<XYChart.Series<Number, Number>> data) {
        this(xAxis, yAxis);
//        xAxis.setAutoRanging(true);
//        this.xAxis.setForceZeroInRange(false);

        setData(data);
    }

    // -------------- METHODS ------------------------------------------------------------------------------------------
    public void setUnifyYAxes(boolean value) {
        unifyYAxes = value;
    }

    /**
     * Called to update and layout the content for the plot
     */
    @Override
    protected void layoutPlotChildren() {
        // we have nothing to layout if no data is present
        addHasData();
        if (getData() == null) {
            return;
        }
        // update candle positions
        int nSeries = getData().size();
        for (int seriesIndex = 0; seriesIndex < getData().size(); seriesIndex++) {
            XYChart.Series<Number, Number> series = getData().get(seriesIndex);
//            System.out.println("series " + series.getName());
            Iterator<XYChart.Data<Number, Number>> iter = getDisplayedDataIterator(series);
            Path seriesPath = null;
            if (series.getNode() instanceof Path) {
                seriesPath = (Path) series.getNode();
                seriesPath.getElements().clear();
            }
            while (iter.hasNext()) {
                XYChart.Data<Number, Number> item = iter.next();
                double yValue = getCurrentDisplayedYValue(item).doubleValue();
                double x = getXAxis().getDisplayPosition(getCurrentDisplayedXValue(item));
                double y = getYAxis().getDisplayPosition(getCurrentDisplayedYValue(item));
                double zeroPos = getYAxis().getZeroPosition();
                Node itemNode = item.getNode();
                ErrorExtraValues extra = (ErrorExtraValues) item.getExtraValue();
                if (itemNode instanceof XYBarWithWhiskers && extra != null) {
                    XYBarWithWhiskers xyBar = (XYBarWithWhiskers) itemNode;

                    double lowPercentile = getYAxis().getDisplayPosition(-extra.getLowPercentile());
                    double highPercentile = getYAxis().getDisplayPosition(-extra.getHighPercentile());
                    // calculate candle width
                    double barWidth = -1;
                    double fullWidth = 1;
                    double stepSize = 1.0;
                    double gap = 0.05;
                    if (getXAxis() instanceof NumberAxis) {
                        NumberAxis xa = (NumberAxis) getXAxis();
                        stepSize = (xa.getDisplayPosition(1.0) - xa.getDisplayPosition(0.0));
                        fullWidth = stepSize / nSeries;
                        barWidth = fullWidth - gap * stepSize / nSeries;
                    }
                    // update candle
//                    System.out.println(x + " " + y + " " + zeroPos + " " + lowPercentile + " " + highPercentile);
                    xyBar.update(zeroPos - y, zeroPos - lowPercentile, zeroPos - highPercentile, barWidth);

//                    xyBar.updateTooltip(item.getYValue().doubleValue(), extra.getClose(), extra.getHigh(), extra.getLow());
                    // position the candle
                    xyBar.setLayoutX(x - fullWidth * (nSeries - 1) / 2 + fullWidth * seriesIndex - barWidth / 2);
                    xyBar.setLayoutX(x - 0.5 * stepSize + gap * stepSize / 2.0 + barWidth * seriesIndex);
                    xyBar.setLayoutY(y);
                    final int seriesIndexFinal = seriesIndex;
                    String seriesName = series.getName();
                    xyBar.setOnMouseClicked(e -> mouseClicked(e, seriesName, seriesIndexFinal, item));

                }
            }
        }
        if (painter != null) {
            painter.addSS(group);
        }

    }

    void mousePressed(MouseEvent e) {
        PyController.mainController.activeChart = this;
    }

    void mouseClickedOnPresenceBar(MouseEvent e, int resNum, double width) {
        boolean appendMode = e.isShiftDown();
        double f = e.getX() / width;
        int seriesIndex = (int) Math.floor(getData().size() * f);
        String seriesName = getData().get(seriesIndex).getName();
        showInfo(seriesName, seriesIndex, resNum, appendMode);
    }

    void mouseClicked(MouseEvent e, String seriesName, int seriesIndex, XYChart.Data<Number, Number> item) {
        boolean appendMode = e.isShiftDown();
        int resNum = (int) Math.round(item.getXValue().doubleValue());
        showInfo(seriesName, seriesIndex, resNum, appendMode);
    }

    void showInfo(String seriesName, int seriesIndex, int resNum, boolean appendMode) {
        String resString = String.valueOf(resNum);
        if (!appendMode) {
            selectedResidues.clear();
            selectedResidues.add(resString);
        } else {
            if (!selectedResidues.contains(resString)) {
                selectedResidues.add(resString);
            }
        }
        currentSeriesName = seriesName;
        String[] seriesNameParts = seriesName.split("\\|");
        String mapName = seriesNameParts[0];
        ResidueProperties resProps = ChartUtil.residueProperties.get(mapName);
        showInfo(resProps, seriesName);
    }

    void showInfo() {
        String[] seriesNameParts = currentSeriesName.split("\\|");
        String mapName = seriesNameParts[0];
        ResidueProperties resProps = ChartUtil.residueProperties.get(mapName);
        showInfo(resProps, currentSeriesName);
    }

    void showInfo(ResidueProperties resProps, String seriesName) {
        Node chartNode = ChartUtil.findNode(this.getScene(), "cpmgchart");
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
        ArrayList<PlotEquation> equations = new ArrayList<>();
        ObservableList<XYChart.Series<Double, Double>> allData = FXCollections.observableArrayList();
        List<ResidueData> resDatas = new ArrayList<>();
        List<int[]> allStates = new ArrayList<>();
        for (ExperimentData expData : resProps.getExperimentData()) {
            for (String residue : residues) {
                resDatas.add(expData.getResidueData(residue));
            }
            String expName = expData.getName();
            if (!ResidueProperties.matchStateString(state, expData.getState())) {
                continue;
            }
            List<XYChart.Series<Double, Double>> data = ChartUtil.getMapData(mapName, expName, residues);
            allData.addAll(data);
            equations.addAll(ChartUtil.getEquations(mapName, residues, equationName, expData.getState(), expData.getField()));
            int[] states = resProps.getStateIndices(0, expData);
            allStates.add(states);
        }
        ((XYChart) chartNode).setData(allData);

        ((PlotData) chartNode).setEquations(equations);
        ((PlotData) chartNode).layoutPlotChildren();
        PyController.mainController.updateTable(resDatas);
        PyController.mainController.updateTableWithPars(mapName, residues, equationName, state, allStates);
    }

    @Override
    protected void dataItemChanged(XYChart.Data<Number, Number> item) {
    }

    @Override
    protected void dataItemAdded(XYChart.Series<Number, Number> series, int itemIndex, XYChart.Data<Number, Number> item) {
        Node candle = createCandle(getData().indexOf(series), item, itemIndex);
        if (shouldAnimate()) {
            candle.setOpacity(0);
            getPlotChildren().add(candle);
            // fade in new candle
            FadeTransition ft = new FadeTransition(Duration.millis(500), candle);
            ft.setToValue(1);
            ft.play();
        } else {
            getPlotChildren().add(candle);
        }
        // always draw average line on top
        if (series.getNode() != null) {
            series.getNode().toFront();
        }
    }

    @Override
    protected void dataItemRemoved(XYChart.Data<Number, Number> item, XYChart.Series<Number, Number> series) {
        final Node candle = item.getNode();
        if (shouldAnimate()) {
            // fade out old candle
            FadeTransition ft = new FadeTransition(Duration.millis(500), candle);
            ft.setToValue(0);
            ft.setOnFinished((ActionEvent actionEvent) -> {
                getPlotChildren().remove(candle);
            });
            ft.play();
        } else {
            getPlotChildren().remove(candle);
        }
    }

    @Override
    protected void seriesAdded(XYChart.Series<Number, Number> series, int seriesIndex) {
        // handle any data already in series
        setTitle(series.getName());
        for (int j = 0; j < series.getData().size(); j++) {
            XYChart.Data item = series.getData().get(j);

            Node candle = createCandle(seriesIndex, item, j);
            if (shouldAnimate()) {
                candle.setOpacity(0);
                getPlotChildren().add(candle);
                // fade in new candle
                FadeTransition ft = new FadeTransition(Duration.millis(500), candle);
                ft.setToValue(1);
                ft.play();
            } else {
                getPlotChildren().add(candle);
            }
        }
        // create series path
//        Path seriesPath = new Path();
//        seriesPath.getStyleClass().setAll("candlestick-average-line", "series" + seriesIndex);
//        series.setNode(seriesPath);
//        getPlotChildren().add(seriesPath);
    }

    @Override
    protected void seriesRemoved(XYChart.Series<Number, Number> series) {
        // remove all candle nodes
        for (XYChart.Data<Number, Number> d : series.getData()) {
            final Node candle = d.getNode();
            if (shouldAnimate()) {
                // fade out old candle
                FadeTransition ft = new FadeTransition(Duration.millis(500), candle);
                ft.setToValue(0);
                ft.setOnFinished((ActionEvent actionEvent) -> {
                    getPlotChildren().remove(candle);
                });
                ft.play();
            } else {
                getPlotChildren().remove(candle);
            }
        }
    }

    /**
     * Create a new Candle node to represent a single data item
     *
     * @param seriesIndex The index of the series the data item is in
     * @param item The data item to create node for
     * @param itemIndex The index of the data item in the series
     * @return New candle node to represent the give data item
     */
    private Node createCandle(int seriesIndex, final XYChart.Data item, int itemIndex) {
        Node xyBar = item.getNode();
        // check if candle has already been created
        if (xyBar instanceof XYBarWithWhiskers) {
            ((XYBarWithWhiskers) xyBar).setSeriesAndDataStyleClasses("series" + seriesIndex, "data" + itemIndex);
        } else {
            xyBar = new XYBarWithWhiskers("series" + seriesIndex, "data" + itemIndex);
            item.setNode(xyBar);
        }
        return xyBar;
    }

    /**
     * This is called when the range has been invalidated and we need to update it. If the axis are auto ranging then we
     * compile a list of all data that the given axis has to plot and call invalidateRange() on the axis passing it that
     * data.
     */
    @Override
    protected void updateAxisRange() {
        // For candle stick chart we need to override this method as we need to let the axis know that they need to be able
        // to cover the whole area occupied by the high to low range not just its center data value
        if (!unifyYAxes) {
            super.updateAxisRange();
            return;
        }
        final Axis<Number> xa = getXAxis();
        final Axis<Number> ya = getYAxis();
        List<Number> xData = null;
        List<Number> yData = null;
        if (xa.isAutoRanging()) {
            xData = new ArrayList<Number>();
        }
        if (ya.isAutoRanging()) {
            yData = new ArrayList<Number>();
        }
        boolean foundPar = false;
        for (String resPropName : ChartUtil.residueProperties.keySet()) {
            ObservableList<XYChart.Series<Double, Double>> data = ChartUtil.getParMapData(resPropName, "best", "0:0:0", "Kex");
            if ((data != null) && (xData != null || yData != null)) {
                foundPar = true;
                for (XYChart.Series<Double, Double> series : data) {
                    for (XYChart.Data<Double, Double> xyData : series.getData()) {
                        if (xData != null) {
                            xData.add(xyData.getXValue());
                        }
                        if (yData != null) {
                            ErrorExtraValues extras = (ErrorExtraValues) xyData.getExtraValue();
                            if (extras != null) {
                                yData.add(xyData.getYValue().doubleValue() + extras.getHighPercentile());
                                yData.add(xyData.getYValue().doubleValue() - extras.getLowPercentile());
                            } else {
                                yData.add(xyData.getYValue());
                            }
                        }
                    }
                }
            }

        }
        if (xData != null) {
            xa.invalidateRange(xData);
        }
        if (yData != null) {
            ya.invalidateRange(yData);
        }
        if (!foundPar) {
            super.updateAxisRange();
        }

    }

    public void setResProps(ResidueProperties resProps) {
        this.resProps = resProps;
    }

    public void addHasData() {
        for (Node node : (ObservableList<Node>) getPlotChildren()) {
            if ((node.getId() != null) && node.getId().equals("presentIndicator")) {
                presenceGroup = (Group) node;
            }
        }
        if (presenceGroup == null) {
            presenceGroup = new Group();
            presenceGroup.setId("presentIndicator");
            getPlotChildren().add(0, presenceGroup);
        }
        if (!presenceGroup.getChildren().isEmpty()) {
            presenceGroup.getChildren().clear();
        }
        if (resProps != null) {

            for (ResidueInfo resInfo : resProps.getResidueMap().values()) {
                if (resInfo == null) {
                    continue;
                }

                int resNum = resInfo.getResNum();
                double x1 = getXAxis().getDisplayPosition(resNum - 0.5) + 1;
                double x2 = getXAxis().getDisplayPosition(resNum + 0.5) - 1;
                double width = x2 - x1;
                double y2 = getYAxis().getHeight();
                Rectangle rect = new Rectangle(width, y2);
                rect.setLayoutX(x1);
                rect.setLayoutY(0.0);
                rect.setFill(Color.LIGHTGRAY);
                rect.setOnMouseClicked(e -> mouseClickedOnPresenceBar(e, resNum, width));

                presenceGroup.getChildren().add(rect);
            }
        }

    }

    public void addSS(ArrayList<SecondaryStructure> secondaryStructures) {
//        painter = new SSPainter(canvas, xAxis, yAxis, secondaryStructures);
        painter = new SSPainter(xAxis, yAxis, secondaryStructures);
        group = null;
        for (Node node : (ObservableList<Node>) getPlotChildren()) {
            if ((node.getId() != null) && node.getId().equals("ssanno")) {
                group = (Group) node;
            }
        }
        if (group == null) {
            group = new Group();
            group.setId("ssanno");
            getPlotChildren().add(1, group);
        }
        layoutPlotChildren();

    }

}

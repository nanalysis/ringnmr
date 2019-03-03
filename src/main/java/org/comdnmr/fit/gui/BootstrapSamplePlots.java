package org.comdnmr.fit.gui;

import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import javafx.beans.Observable;
import javafx.event.ActionEvent;
import javafx.fxml.FXML;
import javafx.scene.Scene;
import javafx.scene.canvas.Canvas;
import javafx.scene.control.Alert;
import javafx.scene.control.Button;
import javafx.scene.control.ChoiceBox;
import javafx.scene.control.Label;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Priority;
import javafx.stage.FileChooser;
import javafx.stage.Stage;
import org.comdnmr.fit.calc.CPMGFitResult;
import org.comdnmr.fit.calc.ResidueInfo;
import org.controlsfx.dialog.ExceptionDialog;
import org.nmrfx.chart.Axis;
import org.nmrfx.chart.DataSeries;
import org.nmrfx.chart.XYCanvasChart;
import org.nmrfx.chart.XYChartPane;
import org.nmrfx.chart.XYValue;
import org.nmrfx.graphicsio.GraphicsIOException;
import org.nmrfx.graphicsio.SVGGraphicsContext;

/**
 *
 * @author brucejohnson
 */
public class BootstrapSamplePlots {

    PyController pyController;
    Stage stage = null;
    XYCanvasChart activeChart;
    BorderPane borderPane = new BorderPane();
    Scene stageScene = new Scene(borderPane, 500, 500);

    ChoiceBox<String> xArrayChoice = new ChoiceBox<>();
    ChoiceBox<String> yArrayChoice = new ChoiceBox<>();

    public BootstrapSamplePlots(PyController controller) {
        pyController = controller;
    }

    public void showMCplot() {
        //Create new Stage for popup window
        if (stage == null) {
            stage = new Stage();
            Label xlabel = new Label("  X Array:  ");
            Label ylabel = new Label("  Y Array:  ");
            //Populate ChoiceBoxes with fitting variable names
            xArrayChoice.getItems().clear();
            yArrayChoice.getItems().clear();
            try {
                xArrayChoice.valueProperty().addListener((Observable x) -> {
                    updateMCplot();
                });
                yArrayChoice.valueProperty().addListener((Observable y) -> {
                    updateMCplot();
                });
            } catch (NullPointerException npEmc1) {
                Alert alert = new Alert(Alert.AlertType.ERROR);
                alert.setContentText("Error: Fit must first be performed.");
                alert.showAndWait();
                return;
            }
            HBox hBox = new HBox();
            HBox.setHgrow(hBox, Priority.ALWAYS);
            Button exportButton = new Button("Export");
            exportButton.setOnAction(e -> exportBarPlotSVGAction(e));

            hBox.getChildren().addAll(exportButton, xlabel, xArrayChoice, ylabel, yArrayChoice);
            //Create the Scatter chart
            XYChartPane chartPane = new XYChartPane();
            activeChart = chartPane.getChart();
            borderPane.setTop(hBox);
            borderPane.setCenter(chartPane);
            stage.setScene(stageScene);
        }
        updateMCPlotChoices();
        stage.show();
        updateMCplot();
    }

    void updateMCplot() {
        Map<String, double[]> simsMap = getMCSimsMap();
        if (simsMap != null) {
            Axis xAxis = activeChart.getXAxis();
            Axis yAxis = activeChart.getYAxis();
            String xElem = xArrayChoice.getValue();
            String yElem = yArrayChoice.getValue();
            if ((xElem != null) && (yElem != null)) {
                xAxis.setLabel(xElem);
                yAxis.setLabel(yElem);
                xAxis.setZeroIncluded(false);
                yAxis.setZeroIncluded(false);
                xAxis.setAutoRanging(true);
                yAxis.setAutoRanging(true);
                DataSeries series = new DataSeries();
                activeChart.getData().clear();
                //Prepare XYChart.Series objects by setting data
                series.getData().clear();
                double[] xValues = simsMap.get(xElem);
                double[] yValues = simsMap.get(yElem);
                if ((xValues != null) && (yValues != null)) {
                    for (int i = 0; i < xValues.length; i++) {
                        series.getData().add(new XYValue(xValues[i], yValues[i]));
                    }
                }
                System.out.println("plot");
                activeChart.getData().add(series);
                activeChart.autoScale(true);
            }
        }
    }

    void updateMCPlotChoices() {
        Map<String, double[]> simsMap = getMCSimsMap();
        xArrayChoice.getItems().clear();
        yArrayChoice.getItems().clear();
        if (simsMap != null) {
            List<String> simParNames = simsMap.keySet().stream().sorted().collect(Collectors.toList());
            xArrayChoice.getItems().addAll(simParNames);
            yArrayChoice.getItems().addAll(simParNames);
            xArrayChoice.setValue(simParNames.get(1));
            yArrayChoice.setValue(simParNames.get(0));
        }
    }

    Map<String, double[]> getMCSimsMap() {
        Map<String, double[]> simsMap = null;
        CPMGFitResult showFitResult = null;
        ResidueInfo currentResInfo = pyController.getResidueInfo();
        String currentEquationName = "";
        if (currentResInfo != null) {
            currentEquationName = pyController.getParametersEquation();
            showFitResult = currentResInfo.getFitResult(currentEquationName);
        }
        if (showFitResult == null) {
            currentEquationName = pyController.simControls.getEquation();
            showFitResult = pyController.getFitResult();
        }
        if (showFitResult != null) {
            simsMap = showFitResult.getSimsMap();
        } else {
            System.out.println("no results");
        }
        stage.setTitle("Monte Carlo Results: " + currentEquationName);
        return simsMap;
    }

    @FXML
    void exportBarPlotSVGAction(ActionEvent event) {
        FileChooser fileChooser = new FileChooser();
        fileChooser.setTitle("Export to SVG");
        fileChooser.setInitialDirectory(pyController.getInitialDirectory());
        File selectedFile = fileChooser.showSaveDialog(null);
        if (selectedFile != null) {
            SVGGraphicsContext svgGC = new SVGGraphicsContext();
            try {
                Canvas canvas = activeChart.getCanvas();
                svgGC.create(true, canvas.getWidth(), canvas.getHeight(), selectedFile.toString());
                exportChart(svgGC);
                svgGC.saveFile();
            } catch (GraphicsIOException ex) {
                ExceptionDialog eDialog = new ExceptionDialog(ex);
                eDialog.showAndWait();
            }
        }
    }

    protected void exportChart(SVGGraphicsContext svgGC) throws GraphicsIOException {
        svgGC.beginPath();
        activeChart.drawChart(svgGC);
    }

}

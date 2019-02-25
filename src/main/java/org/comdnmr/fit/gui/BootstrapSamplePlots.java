package org.comdnmr.fit.gui;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import javafx.beans.Observable;
import javafx.scene.Scene;
import javafx.scene.chart.Axis;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.ScatterChart;
import javafx.scene.chart.XYChart;
import javafx.scene.control.Alert;
import javafx.scene.control.ChoiceBox;
import javafx.scene.control.Label;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Priority;
import javafx.stage.Stage;
import org.comdnmr.fit.calc.CPMGFitResult;
import org.comdnmr.fit.calc.ResidueInfo;

/**
 *
 * @author brucejohnson
 */
public class BootstrapSamplePlots {

    PyController pyController;
    Stage stage = null;
    XYChart activeChart;
    BorderPane borderPane = new BorderPane();
    Scene stageScene = new Scene(borderPane, 500, 500);

    ChoiceBox<String> xArrayChoice = new ChoiceBox<>();
    ChoiceBox<String> yArrayChoice = new ChoiceBox<>();

    public BootstrapSamplePlots(PyController controller) {
        pyController = controller;
    }

    public void showMCplot() {
        //Create new Stage for popup window
        if ((stage == null) || !stage.isShowing()) {
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
            hBox.getChildren().addAll(xlabel, xArrayChoice, ylabel, yArrayChoice);
            //Create the Scatter chart
            NumberAxis MCxAxis = new NumberAxis();
            NumberAxis MCyAxis = new NumberAxis();
            ScatterChart<Number, Number> MCchart = new ScatterChart(MCxAxis, MCyAxis);
            MCxAxis.setAutoRanging(true);
            MCxAxis.setForceZeroInRange(false);
            MCyAxis.setAutoRanging(true);
            MCyAxis.setForceZeroInRange(false);
            activeChart = MCchart;
            activeChart.setLegendVisible(false);
            borderPane.setTop(hBox);
            borderPane.setCenter(activeChart);
            stage.setScene(stageScene);
        }
        updateMCPlotChoices();
        stage.show();
        updateMCplot();
    }

    void updateMCplot() {
        Map<String, double[]> simsMap = getMCSimsMap();
        if (simsMap != null) {
            Axis<Number> xAxis = activeChart.getXAxis();
            Axis<Number> yAxis = activeChart.getYAxis();
            String xElem = xArrayChoice.getValue();
            String yElem = yArrayChoice.getValue();
            xAxis.setLabel(xElem);
            yAxis.setLabel(yElem);
            XYChart.Series<Number, Number> series = new XYChart.Series();
            activeChart.getData().clear();
            activeChart.getData().add(series);
            //Prepare XYChart.Series objects by setting data
            series.getData().clear();
            double[] xValues = simsMap.get(xElem);
            double[] yValues = simsMap.get(yElem);
            if ((xValues != null) && (yValues != null)) {
                for (int i = 0; i < xValues.length; i++) {
                    series.getData().add(new XYChart.Data(xValues[i], yValues[i]));
                }
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

}

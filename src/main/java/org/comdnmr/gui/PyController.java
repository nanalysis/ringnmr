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

import de.jensd.fx.glyphs.GlyphsDude;
import de.jensd.fx.glyphs.fontawesome.FontAwesomeIcon;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import javafx.fxml.FXML;
import javafx.fxml.Initializable;
import java.util.ResourceBundle;
import java.net.URL;
import java.nio.file.FileSystems;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import javafx.application.Platform;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.event.ActionEvent;
import javafx.event.Event;
import javafx.geometry.Bounds;
import javafx.scene.Node;
import javafx.scene.SnapshotParameters;
import javafx.scene.control.CheckBox;
import javafx.scene.control.ChoiceBox;
import javafx.scene.control.Menu;
import javafx.scene.control.MenuItem;
import javafx.scene.control.TabPane;
import javafx.scene.control.TableColumn;
import javafx.scene.control.TableView;
import javafx.scene.control.cell.PropertyValueFactory;
import javafx.scene.image.WritableImage;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.VBox;
import javafx.scene.paint.Color;
import javafx.scene.shape.Circle;
import javafx.stage.FileChooser;
import javafx.stage.Stage;
import javax.imageio.ImageIO;
import javax.script.ScriptException;
import org.comdnmr.eqnfit.CPMGFitter;
import org.comdnmr.eqnfit.FitResult;
import org.comdnmr.eqnfit.CurveFit;
import org.comdnmr.data.DataIO;
import org.comdnmr.eqnfit.EquationFitter;
import org.comdnmr.eqnfit.EquationType;
import org.comdnmr.eqnfit.ExpFitter;
import org.comdnmr.data.Experiment;
import org.comdnmr.data.ExperimentData;
import org.controlsfx.control.PropertySheet;
import org.comdnmr.eqnfit.ParValueInterface;
import org.comdnmr.eqnfit.PlotEquation;
import org.comdnmr.util.ProcessingStatus;
import org.comdnmr.fit.ResidueFitter;
import org.comdnmr.data.ExperimentResult;
import org.comdnmr.data.ExperimentSet;
import org.controlsfx.control.StatusBar;
import org.comdnmr.eqnfit.CESTFitter;
import java.text.DecimalFormat;
import java.util.Collection;
import java.util.Collections;
import java.util.Optional;
import javafx.beans.property.SimpleStringProperty;
import java.util.Random;
import javafx.animation.PauseTransition;
import javafx.embed.swing.SwingFXUtils;
import javafx.fxml.FXMLLoader;
import javafx.print.PrinterJob;
import javafx.scene.Scene;
import javafx.scene.canvas.Canvas;
import javafx.scene.canvas.GraphicsContext;
import javafx.scene.control.Alert;
import javafx.scene.control.SplitPane;
import javafx.scene.control.Label;
import javafx.scene.control.TextField;
import javafx.scene.control.Button;
import javafx.scene.control.ButtonType;
import javafx.scene.control.ContentDisplay;
import javafx.scene.control.ToolBar;
import javafx.scene.input.MouseEvent;
import javafx.scene.layout.Pane;
import javafx.scene.layout.StackPane;
import javafx.stage.StageStyle;
import javafx.util.Duration;
import org.comdnmr.data.DoubleArrayExperiment;
import org.comdnmr.util.CoMDPreferences;
import org.comdnmr.modelfree.CorrelationTime;
import org.comdnmr.eqnfit.FitFunction;
import org.comdnmr.eqnfit.R1RhoFitter;
import static org.comdnmr.gui.MainApp.preferencesController;
import static org.comdnmr.gui.MainApp.console;
import static org.comdnmr.gui.MainApp.primaryStage;
import org.comdnmr.util.CoMDOptions;
import org.comdnmr.utils.NMRFxClient;
import org.controlsfx.dialog.ExceptionDialog;
import org.nmrfx.chart.Axis;
import org.nmrfx.chart.DataSeries;
import org.nmrfx.chart.XYEValue;
import org.nmrfx.chart.XYValue;
import org.nmrfx.chemistry.InvalidMoleculeException;
import org.nmrfx.chemistry.MoleculeBase;
import org.nmrfx.chemistry.MoleculeFactory;
import org.nmrfx.chemistry.relax.RelaxationData.relaxTypes;
import org.nmrfx.chemistry.io.MoleculeIOException;
import org.nmrfx.graphicsio.GraphicsIOException;
import org.nmrfx.graphicsio.SVGGraphicsContext;
import org.nmrfx.peaks.InvalidPeakException;
import org.nmrfx.star.ParseException;
import org.comdnmr.eqnfit.NOEFit;
import org.comdnmr.modelfree.FitModel;
import org.nmrfx.chemistry.relax.RelaxationValues;

public class PyController implements Initializable {

    public static PyController mainController;
    ResidueChart activeChart;
    Stage stage;
    List<ResidueChart> barCharts = new ArrayList<>();

    SSPainter ssPainter = null;
    @FXML
    SSRegion ssregion;

    @FXML
    StackPane stackPane;
    @FXML
    XYPlotDataPane chartPane;

    @FXML
    Button nmrFxPeakButton;
    @FXML
    TableView resInfoTable;
    @FXML
    TableView parameterTable;
    @FXML
    ChoiceBox<String> equationChoice;
    @FXML
    TabPane parTabPane;
    @FXML
    PropertySheet propertySheet;
    @FXML
    Label aicLabel;
    @FXML
    Label rmsLabel;
    @FXML
    Label rChiSqLabel;

    @FXML
    Pane chartBox;
    @FXML
    StatusBar statusBar;
    Circle statusCircle;
    @FXML
    Menu chartMenu;
    @FXML
    Menu axisMenu;
    @FXML
    BorderPane simPane;

    @FXML
    ChoiceBox<String> simChoice;
    @FXML
    CheckBox sliderGuessCheckBox;
    @FXML
    CheckBox calcErrorsCheckBox;
    @FXML
    TextField nSimTextField;

    @FXML
    TextField xLowerBoundTextField;
    @FXML
    TextField xUpperBoundTextField;
    @FXML
    TextField xTickTextField;
    @FXML
    TextField yLowerBoundTextField;
    @FXML
    TextField yUpperBoundTextField;
    @FXML
    TextField yTickTextField;
    @FXML
    CheckBox scalePlot;

    @FXML
    Button setBoundsButton;
    @FXML
    CheckBox autoscaleXCheckBox;
    @FXML
    CheckBox autoscaleYCheckBox;
    @FXML
    CheckBox zeroXCheckBox;
    @FXML
    CheckBox zeroYCheckBox;

    EquationControls simControls;

    @FXML
    TextField genDataNPtsTextField;
    @FXML
    TextField genDataXLBTextField;
    @FXML
    TextField genDataXUBTextField;
    @FXML
    TextField genDataSDevTextField;
    @FXML
    TextField genDataXValTextField;

    @FXML
    SplitPane splitPane;
    @FXML
    ChoiceBox<String> t1Choice;
    @FXML
    ChoiceBox<String> t2Choice;
    @FXML
    TextField r1MedianField;
    @FXML
    TextField r2MedianField;
    @FXML
    TextField tauCalcField;

    BootstrapSamplePlots bootstrapSamplePlots = null;
    InputDataInterface inputDataInterface = null;
    NMRFxClient cl;

    ResidueFitter residueFitter;
    List<String> fittingResidues = new ArrayList<>();
    boolean simulate = true;
    ChartInfo chartInfo = new ChartInfo();

    FitResult fitResult;
    PlotData xychart;
    Canvas barPlotCanvas = new Canvas();

    @FXML
    ToolBar navigatorToolBar = new ToolBar();

    static Random rand = new Random();
    File initialDir = null;

    @FXML
    private void pyAction(ActionEvent event) {
        Node node = (Node) event.getSource();
        String id = node.getId();
        try {
            MainApp.engine.eval("print('howdy " + id + "')");
            MainApp.engine.eval("onAction('" + id + "')");
        } catch (ScriptException sE) {

        }
        //MainApp.interpreter.exec("onAction(" + node + ")");
    }

    @FXML
    public void displayEquation(ActionEvent event) {
        try {
            simControls.simSliderAction("");
        } catch (NullPointerException npE) {
            Alert alert = new Alert(Alert.AlertType.ERROR);
            alert.setContentText("Error: Residue must be selected in display equation");
            alert.showAndWait();
            return;
        }

    }

    public static PyController create(Stage stage) {
        FXMLLoader loader = new FXMLLoader(PyController.class.getResource("/fxml/RINGScene.fxml"));
        PyController controller = null;

        if (stage == null) {
            stage = new Stage(StageStyle.DECORATED);
        }

        try {
            Scene scene = new Scene((Pane) loader.load());
            stage.setScene(scene);
            scene.getStylesheets().add("/styles/Styles.css");

            controller = loader.<PyController>getController();
            controller.stage = stage;
            //controllers.add(controller);
            stage.setTitle("RING NMR Dynamics");
            stage.show();
        } catch (IOException ioE) {

            ioE.printStackTrace();
            System.out.println(ioE.getMessage());
        }
        return controller;
    }

    @Override
    public void initialize(URL url, ResourceBundle rb) {
        mainController = this;
        //propertySheet.setPropertyEditorFactory(new NvFxPropertyEditorFactory());
//        propertySheet.setMode(PropertySheet.Mode.NAME);
//        propertySheet.setModeSwitcherVisible(false);
//        propertySheet.setSearchBoxVisible(false);
        makeAxisMenu();
//        simControls = new CPMGControls();
        if (getFittingMode().equals("cpmg")) {
            simControls = new CPMGControls();
            xLowerBoundTextField.setText("0.0");
            xUpperBoundTextField.setText("1100.0");
            yLowerBoundTextField.setText("0.0");
            yUpperBoundTextField.setText("65.0");
            xTickTextField.setText("100.0");
            yTickTextField.setText("5.0");
            genDataNPtsTextField.setDisable(true);
            genDataXLBTextField.setDisable(true);
            genDataXUBTextField.setDisable(true);
            genDataXValTextField.setDisable(false);
            StringBuilder sBuilder = new StringBuilder();
            for (double xval : getFitter().getSimXDefaults()) {
                sBuilder.append(String.valueOf(xval));
                sBuilder.append(" ");
            }
            genDataXValTextField.setText(sBuilder.toString());
        } else if (getFittingMode().equals("cest")) {
            simControls = new CESTControls();
            xLowerBoundTextField.setText("-20.0");
            xUpperBoundTextField.setText("20.0");
            if (hasExperimentSet()) {
                if (getCurrentExperimentSet().getExperimentData() != null) {
                    DoubleArrayExperiment expData = (DoubleArrayExperiment) getCurrentExperimentSet().getExperimentData().
                            stream().findFirst().get();
                    double[] xVals = expData.getXVals();
                    xLowerBoundTextField.setText(String.valueOf(Math.floor(xVals[1] / 2) * 2));
                    xUpperBoundTextField.setText(String.valueOf(Math.ceil(xVals[xVals.length - 1] / 2) * 2));
                }
            }
            yLowerBoundTextField.setText("0.0");
            yUpperBoundTextField.setText("1.0");
            xTickTextField.setText("1.0");
            yTickTextField.setText("0.25");
            genDataNPtsTextField.setDisable(false);
            genDataNPtsTextField.setText("100");
            genDataXLBTextField.setDisable(false);
            genDataXLBTextField.setText("-8.0");
            genDataXUBTextField.setDisable(false);
            genDataXUBTextField.setText("8.0");
            genDataXValTextField.setDisable(true);
            genDataXValTextField.setText("");
        } else if (getFittingMode().equals("exp")) {
            simControls = new ExpControls();
            xLowerBoundTextField.setText("0.0");
            xUpperBoundTextField.setText("1.25");
            yLowerBoundTextField.setText("0.0");
            yUpperBoundTextField.setText("100.0");
            xTickTextField.setText("0.25");
            yTickTextField.setText("10.0");
            genDataNPtsTextField.setDisable(true);
            genDataXLBTextField.setDisable(true);
            genDataXUBTextField.setDisable(true);
            genDataXValTextField.setDisable(false);
            StringBuilder sBuilder = new StringBuilder();
            for (double xval : getFitter().getSimXDefaults()) {
                sBuilder.append(String.valueOf(xval));
                sBuilder.append(" ");
            }
            genDataXValTextField.setText(sBuilder.toString());
        } else if (getFittingMode().equals("noe")) {
            simControls = new NOEControls();
            xLowerBoundTextField.setText("-0.5");
            xUpperBoundTextField.setText("1.5");
            yLowerBoundTextField.setText("0.0");
            yUpperBoundTextField.setText("1.25");
            xTickTextField.setText("0.25");
            yTickTextField.setText("0.25");
            genDataNPtsTextField.setDisable(true);
            genDataXLBTextField.setDisable(true);
            genDataXUBTextField.setDisable(true);
            genDataXValTextField.setDisable(false);
            StringBuilder sBuilder = new StringBuilder();
            for (double xval : getFitter().getSimXDefaults()) {
                sBuilder.append(String.valueOf(xval));
                sBuilder.append(" ");
            }
            genDataXValTextField.setText(sBuilder.toString());
        } else if (getFittingMode().equals("r1rho")) {
            simControls = new R1RhoControls();
            xLowerBoundTextField.setText("-20.0");
            xUpperBoundTextField.setText("20.0");
            if (hasExperimentSet()) {
                if (getCurrentExperimentSet().getExperimentData() != null) {
                    DoubleArrayExperiment expData = (DoubleArrayExperiment) getCurrentExperimentSet().getExperimentData().
                            stream().findFirst().get();
                    double[] xVals = expData.getXVals();
                    xLowerBoundTextField.setText(String.valueOf(Math.floor(xVals[1] / 2) * 2));
                    xUpperBoundTextField.setText(String.valueOf(Math.ceil(xVals[xVals.length - 1] / 2) * 2));
                }
            }
            yLowerBoundTextField.setText("0.0");
            yUpperBoundTextField.setText("50.0");
            xTickTextField.setText("1.0");
            yTickTextField.setText("0.25");
            genDataNPtsTextField.setDisable(false);
            genDataNPtsTextField.setText("100");
            genDataXLBTextField.setDisable(false);
            genDataXLBTextField.setText("-8.0");
            genDataXUBTextField.setDisable(false);
            genDataXUBTextField.setText("8.0");
            genDataXValTextField.setDisable(true);
            genDataXValTextField.setText("");
        } else {
            System.out.println("Error: no fitting mode selected.");
        }
        VBox vBox = simControls.makeControls(mainController);
        simPane.centerProperty().set(vBox);
        CoMDOptions options = new CoMDOptions(true);
        residueFitter = new ResidueFitter(options, this::updateFitProgress, this::updateStatus);
        statusCircle = new Circle(10);
        statusBar.getLeftItems().add(statusCircle);
        updateEquationChoices(getFittingMode());
        equationChoice.valueProperty().addListener(e -> {
            equationAction();
        });

        simChoice.getItems().add("Simulate CPMG");
        simChoice.getItems().add("Simulate EXP");
        simChoice.getItems().add("Simulate CEST");
        simChoice.getItems().add("Simulate R1Rho");
        simChoice.setValue("Simulate CPMG");
        simChoice.valueProperty().addListener(s -> {
            simAction();
        });

        splitPane.setDividerPositions(0.4, 0.7);

        setBoundsButton.setOnAction(this::setBounds);

        initResidueNavigator();
        calcErrorsCheckBox.selectedProperty().addListener(e -> FitFunction.setCalcError(calcErrorsCheckBox.isSelected()));
        calcErrorsCheckBox.setSelected(true);
        autoscaleXCheckBox.selectedProperty().addListener(e -> autoscaleX(autoscaleXCheckBox.isSelected()));
        autoscaleXCheckBox.setSelected(false);
        autoscaleYCheckBox.selectedProperty().addListener(e -> autoscaleY(autoscaleYCheckBox.isSelected()));
        autoscaleYCheckBox.setSelected(false);
        zeroXCheckBox.selectedProperty().addListener(e -> includeXZero(zeroXCheckBox.isSelected()));
        zeroXCheckBox.setSelected(false);
        zeroYCheckBox.selectedProperty().addListener(e -> includeYZero(zeroYCheckBox.isSelected()));
        zeroYCheckBox.setSelected(false);
        nmrFxPeakButton.setDisable(true);
        nmrFxPeakButton.setOnAction(e -> nmrFxMessage(e));
        xychart = (PlotData) chartPane.getChart();
        chartPane.widthProperty().addListener(e -> resizeXYPlotCanvas());
        chartPane.heightProperty().addListener(e -> resizeXYPlotCanvas());
        chartBox.widthProperty().addListener(e -> resizeBarPlotCanvas());
        chartBox.heightProperty().addListener(e -> resizeBarPlotCanvas());
        chartBox.getChildren().add(barPlotCanvas);
        addChart();
        barPlotCanvas.setOnMouseClicked(e -> mouseClickedOnBarCanvas(e));
//        mainController.setOnHidden(e -> Platform.exit());
        PauseTransition logoTransition = new PauseTransition(Duration.seconds(5));
        logoTransition.setOnFinished(e -> removeLogo());
        logoTransition.play();
        chartMenu.setOnShowing(e -> {
            makeAxisMenu();
        });
    }

    public Stage getStage() {
        return stage;
    }

    void removeLogo() {
        stackPane.getChildren().remove(1);
    }

    void resizeXYPlotCanvas() {
        Canvas canvas = xychart.getCanvas();
        GraphicsContext gC = canvas.getGraphicsContext2D();
        canvas.setWidth(chartPane.getWidth());
        canvas.setHeight(chartPane.getHeight());
        gC.clearRect(0, 0, canvas.getWidth(), canvas.getHeight());
        xychart.setWidth(canvas.getWidth());
        xychart.setHeight(canvas.getHeight());
    }

    void resizeBarPlotCanvas() {
        double width = chartBox.getWidth();
        double height = chartBox.getHeight();
        barPlotCanvas.setWidth(width);
        barPlotCanvas.setHeight(height);
        GraphicsContext gC = barPlotCanvas.getGraphicsContext2D();
        gC.clearRect(0, 0, barPlotCanvas.getWidth(), barPlotCanvas.getHeight());
        if (ssPainter != null) {
            height -= ssPainter.getHeight();
        }
        double chartHeight = height / barCharts.size();
        double yPos = 0.0;
        for (ResidueChart residueChart : barCharts) {
            residueChart.setWidth(width);
            residueChart.setHeight(chartHeight);
            residueChart.setYPos(yPos);
            residueChart.autoScale(false);
            residueChart.drawChart();
            yPos += chartHeight;
        }
        refreshResidueCharts();
    }

    void setResidueChartBorder() {
        double maxLeftBorder = 0.0;
        for (ResidueChart residueChart : barCharts) {
            double[] borders = residueChart.getMinBorders();
            maxLeftBorder = Math.max(borders[0], maxLeftBorder);
        }
        for (ResidueChart residueChart : barCharts) {
            residueChart.setMinLeftBorder(maxLeftBorder);
        }
    }

    void refreshResidueCharts() {
        setResidueChartBorder();
        GraphicsContext gC = barPlotCanvas.getGraphicsContext2D();
        gC.clearRect(0, 0, barPlotCanvas.getWidth(), barPlotCanvas.getHeight());
        for (ResidueChart residueChart : barCharts) {
            residueChart.drawChart();
        }
        if (ssPainter != null) {
            double ssHeight = ssPainter.getHeight();
            Axis axis = activeChart.xAxis;
            ssPainter.paintSS(axis.getXOrigin(), barPlotCanvas.getHeight() - ssHeight,
                    axis.getWidth(),
                    axis.getLowerBound(), axis.getUpperBound());
        }
    }

    void mouseClickedOnBarCanvas(MouseEvent e) {
        double x = e.getX();
        double y = e.getY();
        for (ResidueChart residueChart : barCharts) {
            Axis xAxis = residueChart.xAxis;
            Axis yAxis = residueChart.yAxis;
            if (residueChart.mouseClicked(e)) {
                activeChart = residueChart;
                break;
            } else if ((x > xAxis.getXOrigin()) && (x < xAxis.getXOrigin() + xAxis.getWidth())) {
                if ((y < yAxis.getYOrigin()) && (y > xAxis.getYOrigin() - yAxis.getHeight())) {
                    activeChart = residueChart;
                    break;
                }
            }
        }
        refreshResidueCharts();
    }

    public void setControls() {
        boolean update = false;
        if (getFittingMode().equals("cpmg") && !(simControls instanceof CPMGControls)) {
            simControls = new CPMGControls();
            update = true;
        } else if (getFittingMode().equals("exp") && !(simControls instanceof ExpControls)) {
            simControls = new ExpControls();
            update = true;
        } else if (getFittingMode().equals("noe") && !(simControls instanceof NOEControls)) {
            simControls = new NOEControls();
            update = true;
        } else if (getFittingMode().equals("cest") && !(simControls instanceof CESTControls)) {
            simControls = new CESTControls();
            update = true;
        } else if (getFittingMode().equals("r1rho") && !(simControls instanceof R1RhoControls)) {
            simControls = new R1RhoControls();
            update = true;
        }
        if (update) {
            VBox vBox = simControls.makeControls(mainController);
            simPane.centerProperty().set(vBox);
            updateEquationChoices(getFittingMode());
        }
        if (hasExperimentSet()) {

            String nucleus = getCurrentExperimentSet().getExperimentData().stream().findFirst().get().getNucleusName();
            simControls.setNucleus(nucleus);
            updateXYChartLabels();
        }
    }

    public void setSimControls() {
        boolean update = false;
        if (getSimMode().equals("cpmg") && !(simControls instanceof CPMGControls)) {
            simControls = new CPMGControls();
            update = true;
            genDataNPtsTextField.setDisable(true);
            genDataXLBTextField.setDisable(true);
            genDataXUBTextField.setDisable(true);
            genDataXValTextField.setDisable(false);
            StringBuilder sBuilder = new StringBuilder();
            for (double xval : getFitter().getSimXDefaults()) {
                sBuilder.append(String.valueOf(xval));
                sBuilder.append(" ");
            }
            genDataXValTextField.setText(sBuilder.toString());
        } else if (getSimMode().equals("exp") && !(simControls instanceof ExpControls)) {
            simControls = new ExpControls();
            update = true;
            genDataNPtsTextField.setDisable(true);
            genDataXLBTextField.setDisable(true);
            genDataXUBTextField.setDisable(true);
            genDataXValTextField.setDisable(false);
            StringBuilder sBuilder = new StringBuilder();
            for (double xval : getFitter().getSimXDefaults()) {
                sBuilder.append(String.valueOf(xval));
                sBuilder.append(" ");
            }
            genDataXValTextField.setText(sBuilder.toString());
        } else if (getSimMode().equals("noe") && !(simControls instanceof NOEControls)) {
            simControls = new NOEControls();
            update = true;
            genDataNPtsTextField.setDisable(true);
            genDataXLBTextField.setDisable(true);
            genDataXUBTextField.setDisable(true);
            genDataXValTextField.setDisable(false);
            StringBuilder sBuilder = new StringBuilder();
            for (double xval : getFitter().getSimXDefaults()) {
                sBuilder.append(String.valueOf(xval));
                sBuilder.append(" ");
            }
            genDataXValTextField.setText(sBuilder.toString());
        } else if (getSimMode().equals("cest") && !(simControls instanceof CESTControls)) {
            simControls = new CESTControls();
            ((CESTControls) simControls).updateDeltaLimits();
            update = true;
            genDataNPtsTextField.setDisable(false);
            genDataNPtsTextField.setText("100");
            genDataXLBTextField.setDisable(false);
            genDataXLBTextField.setText("-8.0");
            genDataXUBTextField.setDisable(false);
            genDataXUBTextField.setText("8.0");
            genDataXValTextField.setDisable(true);
            genDataXValTextField.setText("");
        } else if (getSimMode().equals("r1rho") && !(simControls instanceof R1RhoControls)) {
            simControls = new R1RhoControls();
            ((R1RhoControls) simControls).updateDeltaLimits();
            update = true;
            genDataNPtsTextField.setDisable(false);
            genDataNPtsTextField.setText("100");
            genDataXLBTextField.setDisable(false);
            genDataXLBTextField.setText("-8.0");
            genDataXUBTextField.setDisable(false);
            genDataXUBTextField.setText("8.0");
            genDataXValTextField.setDisable(true);
            genDataXValTextField.setText("");
        }
        if (update) {
            updateEquationChoices(getSimMode());
            VBox vBox = simControls.makeControls(mainController);
            simPane.centerProperty().set(vBox);
        }
    }

    void updateEquationChoices() {
        updateEquationChoices(getFittingMode());
    }

    void updateEquationChoices(String mode) {
        if (mode.equals("cpmg")) {
            simControls.updateEquations(equationChoice, CPMGFitter.getEquationNames());
        } else if (mode.equals("exp")) {
            simControls.updateEquations(equationChoice, ExpFitter.getEquationNames());
        } else if (mode.equals("noe")) {
            simControls.updateEquations(equationChoice, NOEFit.getEquationNames());
        } else if (mode.equals("cest")) {
            if (!CESTFitter.getEquationNames().isEmpty()) {
                simControls.updateEquations(equationChoice, CESTFitter.getEquationNames());
            } else {
                Alert alert = new Alert(Alert.AlertType.ERROR, "At least one CEST equation must be selected in Preferences.");
                alert.showAndWait();
            }
        } else if (mode.equals("r1rho")) {
            if (!R1RhoFitter.getEquationNames().isEmpty()) {
                simControls.updateEquations(equationChoice, R1RhoFitter.getEquationNames());
            } else {
                Alert alert = new Alert(Alert.AlertType.ERROR, "At least one R1RHO equation must be selected in Preferences.");
                alert.showAndWait();
            }
        }
    }

    void initResidueNavigator() {

        String iconSize = "12px";
        String fontSize = "7pt";
        ArrayList<Button> buttons = new ArrayList<>();
        Button bButton;

        bButton = GlyphsDude.createIconButton(FontAwesomeIcon.FAST_BACKWARD, "", iconSize, fontSize, ContentDisplay.GRAPHIC_ONLY);
        bButton.setOnAction(e -> firstResidue(e));
        buttons.add(bButton);
        bButton = GlyphsDude.createIconButton(FontAwesomeIcon.BACKWARD, "", iconSize, fontSize, ContentDisplay.GRAPHIC_ONLY);
        bButton.setOnAction(e -> previousResidue(e));
        buttons.add(bButton);
        bButton = GlyphsDude.createIconButton(FontAwesomeIcon.FORWARD, "", iconSize, fontSize, ContentDisplay.GRAPHIC_ONLY);
        bButton.setOnAction(e -> nextResidue(e));
        buttons.add(bButton);
        bButton = GlyphsDude.createIconButton(FontAwesomeIcon.FAST_FORWARD, "", iconSize, fontSize, ContentDisplay.GRAPHIC_ONLY);
        bButton.setOnAction(e -> lastResidue(e));
        buttons.add(bButton);

        navigatorToolBar.getItems().addAll(buttons);
    }

    public void previousResidue(ActionEvent event) {
        List<ExperimentResult> resInfo = getCurrentExperimentSet().getExperimentResults();
        List resNums = new ArrayList<>();
        for (int i = 0; i < resInfo.size(); i++) {
            resNums.add(resInfo.get(i).getResNum());
        }
        Collections.sort(resNums);
        if (chartInfo.hasResidue()) {
            int resIndex = resNums.indexOf(chartInfo.getResNum());
            resIndex--;
            if (resIndex <= 0) {
                resIndex = 0;
            }
            int res = (int) resNums.get(resIndex);
            ResidueChart chart = getActiveChart();
            chart.showInfo(chart.currentSeriesName, 0, res, false);
        }
    }

    public void firstResidue(ActionEvent event) {
        if (chartInfo.hasResidue()) {
            int res = ChartUtil.minRes;
            ResidueChart chart = getActiveChart();
            chart.showInfo(chart.currentSeriesName, 0, res, false);
        }
    }

    public void nextResidue(ActionEvent event) {
        List<ExperimentResult> resInfo = getCurrentExperimentSet().getExperimentResults();
        List resNums = new ArrayList<>();
        for (int i = 0; i < resInfo.size(); i++) {
            resNums.add(resInfo.get(i).getResNum());
        }
        Collections.sort(resNums);
        if (chartInfo.hasResidue()) {
            int resIndex = resNums.indexOf(chartInfo.getResNum());
            resIndex++;
            if (resIndex >= resNums.size()) {
                resIndex = resNums.size() - 1;
            }
            int res = (int) resNums.get(resIndex);
            ResidueChart chart = getActiveChart();
            chart.showInfo(chart.currentSeriesName, 0, res, false);
        }
    }

    public void lastResidue(ActionEvent event) {
        if (chartInfo.hasResidue()) {
            int res = ChartUtil.maxRes;
            ResidueChart chart = getActiveChart();
            chart.showInfo(chart.currentSeriesName, 0, res, false);
        }
    }

    public void loadParameterFile(Event e) {
        FileChooser fileChooser = new FileChooser();
        fileChooser.setTitle("Open YAML File");
        fileChooser.getExtensionFilters().addAll(new FileChooser.ExtensionFilter("Yaml File", "*.yaml", "*.yml"));
        Stage stage = MainApp.primaryStage;
        File file = fileChooser.showOpenDialog(stage);
        if (file != null) {
            if (activeChart != null) {
                clearChart();
                chartInfo.clear();
                simulate = false;
                fitResult = null;
            }
            ChartUtil.loadParameters(file.toString());
        }
        clearSecondaryStructure();
    }

    @FXML
    public void loadPDBFile(Event e) throws MoleculeIOException, ParseException {
        FileChooser fileChooser = new FileChooser();
        fileChooser.setTitle("Open PDB File");
        fileChooser.getExtensionFilters().addAll(new FileChooser.ExtensionFilter("PDB File", "*.pdb"));
        Stage stage = MainApp.primaryStage;
        File file = fileChooser.showOpenDialog(stage);
        String type = "pdb";
        if (file != null) {
            ChartUtil.loadMoleculeFile(file.toString(), type);
        }
    }

    @FXML
    public void loadSTARFile(Event e) throws MoleculeIOException, ParseException {
        FileChooser fileChooser = new FileChooser();
        fileChooser.setTitle("Open STAR File");
        fileChooser.getExtensionFilters().addAll(new FileChooser.ExtensionFilter("STAR File", "*.str"));
        Stage stage = MainApp.primaryStage;
        File file = fileChooser.showOpenDialog(stage);
        String type = "star";
        if (file != null) {
            ChartUtil.loadMoleculeFile(file.toString(), type);
        }
    }

    @FXML
    public void loadNEFFile(Event e) throws MoleculeIOException, ParseException {
        FileChooser fileChooser = new FileChooser();
        fileChooser.setTitle("Open NEF File");
        fileChooser.getExtensionFilters().addAll(new FileChooser.ExtensionFilter("NEF File", "*.nef"));
        Stage stage = MainApp.primaryStage;
        File file = fileChooser.showOpenDialog(stage);
        String type = "nef";
        if (file != null) {
            ChartUtil.loadMoleculeFile(file.toString(), type);
        }
    }

    @FXML
    public void loadCIFFile(Event e) throws MoleculeIOException, ParseException {
        FileChooser fileChooser = new FileChooser();
        fileChooser.setTitle("Open CIF File");
        fileChooser.getExtensionFilters().addAll(new FileChooser.ExtensionFilter("CIF File", "*.cif"));
        Stage stage = MainApp.primaryStage;
        File file = fileChooser.showOpenDialog(stage);
        String type = "cif";
        if (file != null) {
            ChartUtil.loadMoleculeFile(file.toString(), type);
        }
    }

    @FXML
    public void loadSecondaryStructure() {
        List<SecondaryStructure> ssValues = SecondaryStructure.loadFromFile();
        if (!ssValues.isEmpty()) {
            ssPainter = new SSPainter(barPlotCanvas, ssValues);
        }
        resizeBarPlotCanvas();
    }

    @FXML
    public void clearSecondaryStructure() {
        if (ssPainter != null) {
            ssPainter.secondaryStructures.clear();
            ssPainter = null;
            resizeBarPlotCanvas();
        }
    }

    @FXML
    public void inputParameters(ActionEvent event) {
        if (inputDataInterface == null) {
            inputDataInterface = new InputDataInterface(this);
        }
        inputDataInterface.inputParameters();
    }

    @FXML
    public void startServer(ActionEvent event) {
        String tempDir = System.getProperty("java.io.tmpdir");
        String userName = System.getProperty("user.name");
        Path path = FileSystems.getDefault().getPath(tempDir, "NMRFx_" + userName + "_port.txt");
        File f = path.toFile(); //new File(tempDir + "/NMRFx_" + userName + "_port.txt");
        int port = 8021;
        if (!f.exists() && !f.isDirectory()) {
            Alert alert = new Alert(Alert.AlertType.WARNING, "NMRFx Server port file " + "\n" + f + "\ndoesn't exist. \nCreate NMRFx server first in NMRFxProcessor GUI.");
            alert.showAndWait();
        } else {
            try {
                BufferedReader reader = new BufferedReader(new FileReader(f));
                String text = null;
                while ((text = reader.readLine()) != null) {
                    port = Integer.parseInt(text);
                }
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            } catch (IOException e) {
                e.printStackTrace();
            }
            cl = NMRFxClient.makeClient(port);
            nmrFxPeakButton.setDisable(false);
        }
    }

    public NMRFxClient getClient() {
        return cl;
    }

    @FXML
    void nmrFxMessage(ActionEvent e) {
        String peakNum = getPeakNumFromTable();
        NMRFxClient cl = PyController.mainController.getClient();
        try {
            String[] peakSplit = peakNum.split("\\.");
            String peakName = peakSplit[0];
            String peakNumber = peakSplit[1];
            String peakString = peakNum;
            if (!peakName.equals("")) {
                peakString = peakNumber + "/" + peakName;
            } else if (peakName.equals("")) {
                peakString = peakNumber;
            }
            cl.sendMessage("showpeak/" + peakString);
        } catch (IOException ioE) {
            System.out.println(ioE.getMessage());
        }
    }

    public void updateXYChartLabels() {
        if ((simControls instanceof CPMGControls)) {
            xychart.setNames("CPMG", "\u03BD (cpmg)", "R2 (\u03BD)", "0");
            xychart.setBounds(0.0, 1100.0, 0.0, 65.0, 100.0, 5.0);
            xLowerBoundTextField.setText("0.0");
            xUpperBoundTextField.setText("1100.0");
            yLowerBoundTextField.setText("0.0");
            yUpperBoundTextField.setText("65.0");
            xTickTextField.setText("100.0");
            yTickTextField.setText("5.0");
        } else if ((simControls instanceof ExpControls)) {
            xychart.setNames("Exp", "Time (s)", "Intensity", "0");
            xychart.setBounds(0.0, 1.25, 0.0, 100, 0.25, 10.0);
            xLowerBoundTextField.setText("0.0");
            xUpperBoundTextField.setText("1.25");
            yLowerBoundTextField.setText("0.0");
            yUpperBoundTextField.setText("100.0");
            xTickTextField.setText("0.25");
            yTickTextField.setText("10.0");
        } else if ((simControls instanceof NOEControls)) {
            xychart.setNames("NOE", "On/Off", "Intensity", "0");
            xychart.setBounds(-0.25, 1.25, 0.0, 1.25, 0.25, 10.0);
            xLowerBoundTextField.setText("00.25");
            xUpperBoundTextField.setText("1.25");
            yLowerBoundTextField.setText("0.0");
            yUpperBoundTextField.setText("1.25");
            xTickTextField.setText("0.25");
            yTickTextField.setText("0.25");
        } else if ((simControls instanceof CESTControls)) {
            xychart.setNames("CEST", "Offset (PPM)", "I(t)/I(0)", "20");
            xychart.setBounds(-20, 20, 0.0, 1.0, 2.0, 0.25);
            xLowerBoundTextField.setText("-20.0");
            xUpperBoundTextField.setText("20.0");
            if (hasExperimentSet()) {
                DoubleArrayExperiment expData = (DoubleArrayExperiment) getCurrentExperimentSet().getExperimentData().
                        stream().findFirst().get();
                double[] xVals = expData.getXVals();
                if (xVals != null) {
                    xychart.setBounds(Math.floor(xVals[1] / 2) * 2, Math.ceil(xVals[xVals.length - 1] / 2) * 2, 0.0, 1.0, 1.0, 0.25);
                    xLowerBoundTextField.setText(String.valueOf(Math.floor(xVals[1] / 2) * 2));
                    xUpperBoundTextField.setText(String.valueOf(Math.ceil(xVals[xVals.length - 1] / 2) * 2));
                }
            }
            yLowerBoundTextField.setText("0.0");
            yUpperBoundTextField.setText("1.0");
            xTickTextField.setText("2.0");
            yTickTextField.setText("0.25");
        } else if ((simControls instanceof R1RhoControls)) {
            xychart.setNames("R1Rho", "Offset (PPM)", "R1Rho", "20");
            xychart.setBounds(-20, 20, 0.0, 50.0, 2.0, 5.0);
            xLowerBoundTextField.setText("-20.0");
            xUpperBoundTextField.setText("20.0");
            if (hasExperimentSet()) {
                DoubleArrayExperiment expData = (DoubleArrayExperiment) getCurrentExperimentSet().getExperimentData().
                        stream().findFirst().get();
                double[] xVals = expData.getXVals();
                if (xVals != null) {
                    xychart.setBounds(Math.floor(xVals[1] / 2) * 2, Math.ceil(xVals[xVals.length - 1] / 2) * 2, 0.0, 1.0, 1.0, 0.25);
                    xLowerBoundTextField.setText(String.valueOf(Math.floor(xVals[1] / 2) * 2));
                    xUpperBoundTextField.setText(String.valueOf(Math.ceil(xVals[xVals.length - 1] / 2) * 2));
                }
            }
            yLowerBoundTextField.setText("0.0");
            yUpperBoundTextField.setText("50.0");
            xTickTextField.setText("2.0");
            yTickTextField.setText("5.0");
        }
    }

    public void setBounds(ActionEvent event) {
        setBounds();
    }

    public void setBounds() {
        try {
            double xLB = Double.parseDouble(xLowerBoundTextField.getText());
            double xUB = Double.parseDouble(xUpperBoundTextField.getText());
            double xTick = Double.parseDouble(xTickTextField.getText());
            double yLB = Double.parseDouble(yLowerBoundTextField.getText());
            double yUB = Double.parseDouble(yUpperBoundTextField.getText());
            double yTick = Double.parseDouble(yTickTextField.getText());
            if (xLB < xUB & yLB < yUB) {
                xychart.setBounds(xLB, xUB, yLB, yUB, xTick, yTick);
            } else {
                Alert alert = new Alert(Alert.AlertType.ERROR);
                alert.setContentText("Error: Upper Bound must be > Lower Bound.");
                alert.showAndWait();
                return;
            }

        } catch (NumberFormatException nfEsb) {
            Alert alert = new Alert(Alert.AlertType.ERROR);
            alert.setContentText("Error: Bound and Interval values must be provided.");
            alert.showAndWait();
            return;
        }

    }

    public void autoscaleBounds(ActionEvent event) {
        double[] bounds = xychart.autoScale(true);
        if (bounds != null) {
            xLowerBoundTextField.setText(Double.toString(bounds[0]));
            xUpperBoundTextField.setText(Double.toString(bounds[1]));
            yLowerBoundTextField.setText(Double.toString(bounds[2]));
            yUpperBoundTextField.setText(Double.toString(bounds[3]));
            if (simControls instanceof CESTControls) {
                ((CESTControls) simControls).updateDeltaLimits(bounds[0], bounds[1]);
            }
        }

    }

    public void autoscaleX(boolean autoX) {
        xychart.xAxis.setAutoRanging(autoX);
    }

    public void autoscaleY(boolean autoY) {
        xychart.yAxis.setAutoRanging(autoY);
    }

    public void includeXZero(boolean xZero) {
        xychart.xAxis.setZeroIncluded(xZero);
    }

    public void includeYZero(boolean yZero) {
        xychart.yAxis.setZeroIncluded(yZero);
    }

    public void updateChartEquations(String equationName, double[] pars, double[] errs, double[] fields) {
        List<GUIPlotEquation> equations = new ArrayList<>();
        String expType = getFittingMode();
        for (int i = 0; i < fields.length; i++) {
            double[] extras = {fields[i] / fields[0]};
            //System.out.println("updateChartEquations got called with extras length = "+extras.length);
            GUIPlotEquation plotEquation = new GUIPlotEquation(expType, equationName, pars, errs, extras);
            equations.add(plotEquation);
        }
        showEquations(equations);
    }

    public void showEquations(List<GUIPlotEquation> equations) {
        xychart.setEquations(equations);
//        Optional<Double> rms = rms();

    }

    public void calcModel1() {
        FitModel fitModel = new FitModel();
        fitModel.testIsoModel();

    }

    public void estimateCorrelationTime() {
        String r1SetName = t1Choice.getValue();
        String r2SetName = t2Choice.getValue();
        ExperimentSet r1Set = ChartUtil.getResidueProperty(r1SetName);
        ExperimentSet r2Set = ChartUtil.getResidueProperty(r2SetName);

        Map<String, Double> result = CorrelationTime.estimateTau(r1Set, r2Set);
        if (!result.isEmpty()) {
            r1MedianField.setText(String.format("%.3f 1/s", result.get("R1")));
            r2MedianField.setText(String.format("%.3f 1/s", result.get("R2")));
            tauCalcField.setText(String.format("%.2f ns", result.get("tau")));
        }
    }

    public void model1Order() {
//        String r1SetName = t1Choice.getValue();
//        String r2SetName = t2Choice.getValue();
//        String noeSetName = noeChoice.getValue();
//        CorrelationTime.fitS(ChartUtil.residueProperties,
//                r1SetName, r2SetName, noeSetName);
    }

    public ResidueChart getActiveChart() {
        return activeChart;
    }

    public ResidueChart addChart() {
        ResidueChart newChart = ResidueChart.buildChart(barPlotCanvas);
        barCharts.add(newChart);
        activeChart = newChart;
        resizeBarPlotCanvas();
        return activeChart;
    }

    public void removeChart(Event e) {
        if ((activeChart != null) && (barCharts.size() > 1)) {
            barCharts.remove(activeChart);
            activeChart = barCharts.get(0);
            resizeBarPlotCanvas();
        }
    }

    public void updateTable(List<ExperimentData> experimentalDataSets) {
        ObservableList<ExperimentData.DataValue> data = FXCollections.observableArrayList();
        for (ExperimentData experimentalData : experimentalDataSets) {
            data.addAll(experimentalData.getDataValues());
        }
        resInfoTable.itemsProperty().setValue(data);

        TableColumn<ExperimentData.DataValue, String> nameColumn = new TableColumn<>("Name");
        TableColumn<ExperimentData.DataValue, String> resColumn = new TableColumn<>("Residue");
        TableColumn<ExperimentData.DataValue, String> resNameColumn = new TableColumn<>("ResName");
        TableColumn<ExperimentData.DataValue, String> errColumn = new TableColumn<>("Error");
        TableColumn<ExperimentData.DataValue, String> peakColumn = new TableColumn<>("Peak");

        nameColumn.setCellValueFactory(new PropertyValueFactory<>("Name"));
        resColumn.setCellValueFactory(new PropertyValueFactory<>("Residue"));
        resNameColumn.setCellValueFactory(new PropertyValueFactory<>("ResName"));
        errColumn.setCellValueFactory(new PropertyValueFactory<>("Error"));
        peakColumn.setCellValueFactory(new PropertyValueFactory<>("Peak"));

        resInfoTable.getColumns().clear();
        resInfoTable.getColumns().addAll(nameColumn, resColumn, resNameColumn, errColumn, peakColumn);

        if (getFittingMode().equals("cpmg")) {
            TableColumn<ExperimentData.DataValue, Double> xColumn = new TableColumn<>("Vcpmg");
            TableColumn<ExperimentData.DataValue, Double> yColumn = new TableColumn<>("Reff");

            xColumn.setCellValueFactory(new PropertyValueFactory<>("X0"));
            yColumn.setCellValueFactory(new PropertyValueFactory<>("Y"));

            resInfoTable.getColumns().clear();
            resInfoTable.getColumns().addAll(nameColumn, resColumn, resNameColumn, xColumn, yColumn, errColumn, peakColumn);
        } else if (getFittingMode().equals("exp")) {
            TableColumn<ExperimentData.DataValue, Double> xColumn = new TableColumn<>("Delay");
            TableColumn<ExperimentData.DataValue, Double> yColumn = new TableColumn<>("Intensity");
//            TableColumn<ExperimentData.DataValue, Double> t1Column = new TableColumn<>("T1");
//            TableColumn<ExperimentData.DataValue, Double> t2Column = new TableColumn<>("T2");
//            TableColumn<ExperimentData.DataValue, Double> t1RhoColumn = new TableColumn<>("T1Rho");

            xColumn.setCellValueFactory(new PropertyValueFactory<>("X0"));
            yColumn.setCellValueFactory(new PropertyValueFactory<>("Y"));
//            t1Column.setCellValueFactory(new PropertyValueFactory<>("T1"));
//            t2Column.setCellValueFactory(new PropertyValueFactory<>("T2"));
//            t1RhoColumn.setCellValueFactory(new PropertyValueFactory<>("T1Rho"));

            resInfoTable.getColumns().clear();
            resInfoTable.getColumns().addAll(nameColumn, resColumn, resNameColumn,
                    //                    t1Column, t2Column, t1RhoColumn,
                    xColumn, yColumn, errColumn, peakColumn);
        } else if (getFittingMode().equals("cest")) {
            TableColumn<ExperimentData.DataValue, Double> x0Column = new TableColumn<>("Offset");
            TableColumn<ExperimentData.DataValue, Double> x1Column = new TableColumn<>("B1 Field");
            TableColumn<ExperimentData.DataValue, Double> yColumn = new TableColumn<>("Intensity");

            x0Column.setCellValueFactory(new PropertyValueFactory<>("X0"));
            x1Column.setCellValueFactory(new PropertyValueFactory<>("X1"));
            yColumn.setCellValueFactory(new PropertyValueFactory<>("Y"));

            resInfoTable.getColumns().clear();
            resInfoTable.getColumns().addAll(nameColumn, resColumn, resNameColumn, x0Column, x1Column, yColumn, errColumn, peakColumn);
        } else if (getFittingMode().equals("r1rho")) {
            TableColumn<ExperimentData.DataValue, Double> x0Column = new TableColumn<>("Offset");
            TableColumn<ExperimentData.DataValue, Double> x1Column = new TableColumn<>("B1 Field");
            TableColumn<ExperimentData.DataValue, Double> yColumn = new TableColumn<>("Intensity");

            x0Column.setCellValueFactory(new PropertyValueFactory<>("X0"));
            x1Column.setCellValueFactory(new PropertyValueFactory<>("X1"));
            yColumn.setCellValueFactory(new PropertyValueFactory<>("Y"));

            resInfoTable.getColumns().clear();
            resInfoTable.getColumns().addAll(nameColumn, resColumn, resNameColumn, x0Column, x1Column, yColumn, errColumn, peakColumn);
        }
    }

//    public void updateTableWithPars(String mapName, String[] residues, String equationName, String state, List<int[]> allStates) {
//        updateTableWithPars(mapName, residues, equationName, state, allStates, true);
//    }
//                updateTableWithPars(currentMapName, currentResidues, equationName, currentState, useStates, false);
    public void updateTableWithPars(ChartInfo chartInfo, boolean savePars) {
        List<ParValueInterface> allParValues = new ArrayList<>();
        if (chartInfo.hasResidues()) {
            for (String resNum : chartInfo.getResidues()) {
                ExperimentResult resInfo = ChartUtil.getResInfo(chartInfo.mapName, String.valueOf(resNum));
                if (resInfo != null) {
                    chartInfo.experimentalResult = resInfo;
                    final String useEquationName;
                    if (chartInfo.equationName.equals("best")) {
                        useEquationName = resInfo.getBestEquationName();
                    } else {
                        useEquationName = chartInfo.equationName;
                    }
                    List<ParValueInterface> parValues = resInfo.getParValues(useEquationName, chartInfo.state);
                    if (resNum.equals(chartInfo.currentResidues[0])) {
                        simControls.updateStates(chartInfo.currentStates);
                        simControls.updateSliders(parValues, useEquationName);
                    }

                    allParValues.addAll(parValues);
                    CurveFit curveSet = chartInfo.experimentalResult.getCurveSet(useEquationName, chartInfo.state.replace("*", "0"));
//                    System.out.println("curv " + useEquationName + " " + state + " " + curveSet);
                    try {
                        String aic = String.format("%.2f", curveSet.getParMap().get("AIC"));
                        String rms = String.format("%.3f", curveSet.getParMap().get("RMS"));
                        String rChiSq = String.format("%.2f", curveSet.getParMap().get("rChiSq"));
                        aicLabel.setText(aic);
                        rmsLabel.setText(rms);
                        rChiSqLabel.setText(rChiSq);
                    } catch (NullPointerException npEaic) {

                    }
                }
            }
            if (savePars) {
//                chartInfo.mapName = mapName;
//                chartInfo.currentResidues = new String[residues.length];
//                System.arraycopy(residues, 0, chartInfo.currentResidues, 0, residues.length);
//                chartInfo.state = state;
//                chartInfo.equationName = equationName;
//                chartInfo.currentStates.clear();
//                chartInfo.currentStates.addAll(allStates);
            }

            updateTableWithPars(allParValues);
        }
    }

    public void updateEquation(String mapName, String[] residues, String equationName) {
        String setEquation = "";
        for (String resNum : residues) {
            final String useEquationName;
            ExperimentResult resInfo = ChartUtil.getResInfo(mapName, resNum);
            if (resInfo != null) {
                if (equationName.equals("best")) {
                    useEquationName = resInfo.getBestEquationName();
                } else {
                    useEquationName = equationName;
                }

                if (setEquation.equals("")) {
                    setEquation = useEquationName;
                } else if (!setEquation.equals(useEquationName)) {
                    setEquation = "+";
                }
            }
        }
        equationChoice.setUserData("disabled");
        equationChoice.setValue(setEquation);
        equationChoice.setUserData(null);
    }

    public void updateTableWithPars(List<ParValueInterface> parValues) {
        DecimalFormat df = new DecimalFormat();
        TableColumn<ParValueInterface, String> residueColumn = new TableColumn<>("Residue");
        TableColumn<ParValueInterface, String> resNameColumn = new TableColumn<>("ResName");
        TableColumn<ParValueInterface, String> stateColumn = new TableColumn<>("State");
        TableColumn<ParValueInterface, String> nameColumn = new TableColumn<>("Name");
        TableColumn<ParValueInterface, String> valueColumn = new TableColumn<>("Value");
        TableColumn<ParValueInterface, String> errorColumn = new TableColumn<>("Error");

        residueColumn.setCellValueFactory(new PropertyValueFactory<>("Residue"));
        resNameColumn.setCellValueFactory(new PropertyValueFactory<>("ResName"));
        stateColumn.setCellValueFactory(new PropertyValueFactory<>("State"));
        nameColumn.setCellValueFactory(new PropertyValueFactory<>("Name"));
        valueColumn.setCellValueFactory(c -> {
            SimpleStringProperty property = new SimpleStringProperty();
            double errValue = c.getValue().getError();
            int nSig = (int) Math.floor(Math.log10(errValue)) - 1;
            nSig = -nSig;
            if (nSig < 0) {
                nSig = 0;
            }
            df.setMaximumFractionDigits(nSig);
            property.setValue(df.format(c.getValue().getValue()));
            return property;
        });
        errorColumn.setCellValueFactory(c -> {
            SimpleStringProperty property = new SimpleStringProperty();
            double errValue = c.getValue().getError();
            int nSig = (int) Math.floor(Math.log10(errValue));
            nSig = -nSig;
            if (nSig < 0) {
                nSig = 0;
            }
            df.setMaximumFractionDigits(nSig);
            property.setValue(df.format(errValue));
            return property;
        });

        parameterTable.getColumns().clear();
        parameterTable.getColumns().addAll(residueColumn, resNameColumn, stateColumn, nameColumn, valueColumn, errorColumn);
        ObservableList<ParValueInterface> data = FXCollections.observableArrayList();

        for (ParValueInterface parValue : parValues) {
            String state = parValue.getState();
            String parName = parValue.getName();
            boolean ok = true;
            // remove redundant parameters 
            if ((state != null) && (parName != null)) {
                if (!state.equals("0:0:0") && (parName.equals("Kex") || parName.equals("dPPM")
                        || parName.equals("pA") || parName.equals("dPPMmin"))) {
                    ok = false;
                }
            }
            if (ok) {
                data.add(parValue);
            }
        }
        parameterTable.itemsProperty().setValue(data);
    }

    public void selectTableRow(String seriesName, int index) {
        parTabPane.getSelectionModel().select(1);
        List<ExperimentData.DataValue> data = resInfoTable.getItems();
        int iRow = 0;
        for (ExperimentData.DataValue dValue : data) {
            if ((dValue.getIndex() == index) && ((dValue.getName() + ":" + dValue.getResidue()).equals(seriesName))) {
                resInfoTable.getSelectionModel().clearAndSelect(iRow);
                resInfoTable.scrollTo(iRow);
                return;
            }
            iRow++;

        }

    }

    public String getPeakNumFromTable() { //getPeakNumFromTable(String seriesName, int index)
        parTabPane.getSelectionModel().select(1);
        List<ExperimentData.DataValue> data = resInfoTable.getItems();
        int iRow = 0;
        int peakNum = 0;
        String name = "";
//        for (ExperimentData.DataValue dValue : data) {
//            if ((dValue.getIndex() == index) && ((dValue.getName() + ":" + dValue.getResidue()).equals(seriesName))) {
//                resInfoTable.getSelectionModel().clearAndSelect(iRow);
//                resInfoTable.scrollTo(iRow);
//            }
//            name = data.get(iRow).getName();
//            peakNum = data.get(iRow).getPeak();
//            iRow++;
//
//        }
        name = data.get(iRow).getName();
        peakNum = data.get(iRow).getPeak();
//        System.out.println("selected row peak num = " + peakNum);
        return name + "." + String.valueOf(peakNum);
    }

    public void clearChart(Event e) {
        clearChart();
    }

    public void clearChart() {
        activeChart.getData().clear();
    }

    public void showRelaxationValues(String setName, String valueName, String parName) {
        List<RelaxationValues> values = ChartUtil.getMolRelaxationValues(setName);
        ObservableList<DataSeries> data = ChartUtil.getRelaxationDataSeries(values, valueName, setName, parName);
        String yLabel = valueName.equalsIgnoreCase(parName) ? parName
                : valueName.toUpperCase() + ": " + parName;

        addSeries(data, setName, yLabel);

    }

    public void setYAxisType(String expMode, String setName, String eqnName, String state, String typeName) {
        ObservableList<DataSeries> data = ChartUtil.getParMapData(setName, eqnName, state, typeName);
        String yLabel = expMode.equalsIgnoreCase(typeName) ? typeName
                : expMode.toUpperCase() + ": " + typeName;
        addSeries(data, setName, yLabel);
    }

    public void addSeries(ObservableList<DataSeries> data, String setName, String yLabel) {
        ResidueChart chart = activeChart;

        chart.getData().addAll(data);
        int i = 0;
        for (DataSeries series : chart.getData()) {
            series.setFill(PlotData.colors[i++ % PlotData.colors.length]);
        }
        if (chart.currentSeriesName.equals("")) {
            chart.currentSeriesName = data.get(0).getName();
        }
        Axis xAxis = chart.xAxis;
        Axis yAxis = chart.yAxis;
        xAxis.setAutoRanging(false);
        yAxis.setAutoRanging(true);
        chart.autoScale(false);
        double xMin = Math.floor((ChartUtil.minRes - 2) / 5.0) * 5.0;
        double xMax = Math.ceil((ChartUtil.maxRes + 2) / 5.0) * 5.0;
        xAxis.setLowerBound(xMin);
        xAxis.setUpperBound(xMax);
//        xAxis.setTickUnit(10);
//        xAxis.setMinorTickCount(5);
//        yAxis.setMinWidth(75.0);
//        yAxis.setPrefWidth(75.0);
        chart.yAxis.setLabel(yLabel);

        chart.setTitle(setName);
        if (chart.getData().size() > 1) {
//            chart.setLegendVisible(true); // fixme
//            chart.setLegendSide(Side.TOP);
        } else {
            // chart.setLegendVisible(false); //fixme
        }
        setCurrentExperimentSet(ChartUtil.getResidueProperty(setName));
        chart.setResProps(getCurrentExperimentSet());
        refreshResidueCharts();
    }

    void makeT1T2Menu() {
        t1Choice.getItems().clear();
        t2Choice.getItems().clear();
        // noeChoice.getItems().clear();
        t1Choice.getItems().add("");
        t2Choice.getItems().add("");
        //noeChoice.getItems().add("");

        Collection<String> setNames = ChartUtil.getMolRelaxationNames();
        for (var setName : setNames) {
            t1Choice.getItems().add(setName);
            t2Choice.getItems().add(setName);
            //  noeChoice.getItems().add(setName);
        }
        setNames = ChartUtil.getResiduePropertyNames();
        for (var setName : setNames) {
            t1Choice.getItems().add(setName);
            t2Choice.getItems().add(setName);
            //  noeChoice.getItems().add(setName);
        }
    }

    void makeAxisMenu() {
        makeT1T2Menu();
        axisMenu.getItems().clear();
        addMoleculeDataToAxisMenu();
        addResiduePropertiesToAxisMenu();
    }

    void addMoleculeDataToAxisMenu() {
        var molResProps = DataIO.getDataFromMolecule();
        ChartUtil.clearMolRelaxationValues();
        for (var entry : molResProps.entrySet()) {
            String setName = entry.getKey();
            ChartUtil.addMolRelaxationValues(setName, molResProps.get(setName));

            Menu cascade = new Menu(setName);
            var values = entry.getValue();
            if (!values.isEmpty()) {
                RelaxationValues value = values.get(0);
                axisMenu.getItems().add(cascade);
                String[] parNames = value.getParNames();
                for (var parName : parNames) {
                    MenuItem cmItem1 = new MenuItem(parName);
                    cmItem1.setOnAction(e -> showRelaxationValues(setName, value.getName(), parName));
                    cascade.getItems().add(cmItem1);
                }

            }

        }
    }

    void addResiduePropertiesToAxisMenu() {
        Collection<String> setNames = ChartUtil.getResiduePropertyNames();
        for (var setName : setNames) {
            var residueProp = ChartUtil.getResidueProperty(setName);
            Menu cascade = new Menu(setName);
            axisMenu.getItems().add(cascade);
            String expMode = residueProp.getExpMode();
            String[] parTypes = getParTypes(residueProp.getExpMode());
            for (String parType : parTypes) {
                Menu cascade2 = new Menu(parType);
                cascade.getItems().add(cascade2);
                ArrayList<String> equationNames = new ArrayList<>();
                equationNames.add("best");
                equationNames.addAll(residueProp.getEquationNames());
                List<String> stateStrings = residueProp.getStateStrings();

                for (String equationName : equationNames) {
                    if ((stateStrings.size() < 2) || parType.equals("RMS") || parType.equals("AIC") || parType.equals("Equation")) {
                        MenuItem cmItem1 = new MenuItem(equationName);
                        cmItem1.setOnAction(e -> setYAxisType(expMode, setName, equationName, "0:0:0", parType));
                        cascade2.getItems().add(cmItem1);
                    } else {
                        boolean isValidPar = false;
                        if (equationName.equals("best")) {
                            isValidPar = true;
                        } else {
                            String[] validPars = getParNames(equationName);
                            for (String validPar : validPars) {
                                if (validPar.equals(parType)) {
                                    isValidPar = true;
                                    break;
                                }
                            }
                        }
                        if (isValidPar) {
                            Menu cascade3 = new Menu(equationName);
                            cascade2.getItems().add(cascade3);
                            residueProp.getStateStrings().stream().forEach(state -> {
                                MenuItem cmItem1 = new MenuItem(state);
                                cmItem1.setOnAction(e -> setYAxisType(expMode, setName, equationName, state, parType));
                                cascade3.getItems().add(cmItem1);

                            });
                        }
                    }

                }
            }
        }
    }

    @FXML
    public void guesses(ActionEvent event) {
        if (!hasExperimentSet()) {
            guessSimData();
        } else {
            try {
                EquationFitter equationFitter = getFitter();
                String[] resNums = {chartInfo.getResNum()};
                equationFitter.setData(getCurrentExperimentSet(), resNums);
                String equationName = simControls.getEquation();
//        System.out.println("guesses eqnFitter = " + equationFitter);
//        System.out.println("guesses resNums = " + resNums);
//        System.out.println("guesses eqnName = " + equationName);
                List<ParValueInterface> guesses = equationFitter.guessPars(equationName);
                if (guesses != null) {
                    simControls.updateSliders(guesses, equationName);
                    simControls.simSliderAction("");
                }
            } catch (NullPointerException npE1) {
                Alert alert = new Alert(Alert.AlertType.ERROR);
                alert.setContentText("Error: Residue must be selected");
                alert.showAndWait();
                return;
            }
        }
    }

    public Optional<Double> rms() {
        Optional<Double> rms = Optional.empty();
        if (hasExperimentSet()) {
            try {
                EquationFitter equationFitter = getFitter();
                if (chartInfo.hasResult()) {
                    String[] resNums = {String.valueOf(chartInfo.getResNum())};
                    equationFitter.setData(getCurrentExperimentSet(), resNums);
                    String equationName = simControls.getEquation();
                    equationFitter.setupFit(equationName);
                    int[][] map = equationFitter.getFitModel().getMap();
                    double[] sliderGuesses = simControls.sliderGuess(equationName, map);
                    rms = Optional.of(equationFitter.rms(sliderGuesses));
                } else {
                    rms = Optional.empty();
                }
            } catch (NullPointerException npE2) {
                rms = Optional.empty();
            }
        }
        return rms;
    }

    @FXML
    public void fitEquation(ActionEvent event) {
//        EquationFitter equationFitter = new CPMGFit();
        fitResult = null;
        try {
            EquationFitter equationFitter = getFitter();
            if (!hasExperimentSet()) {
                fitSimData();
            } else {
                String[] resNums = {chartInfo.getResNum()};
                equationFitter.setData(getCurrentExperimentSet(), resNums);
                String equationName = simControls.getEquation();
                equationFitter.setupFit(equationName);
                int[][] map = equationFitter.getFitModel().getMap();
//            for (int i=0; i<map.length; i++) {
//                for (int j=0; j<map[i].length; j++) {
//                    System.out.println("map " + i + " " + j + " " + map[i][j]);
//                }
//            }
//            System.out.println("getFitModel = " + equationFitter.getFitModel());
//            System.out.println("fitEqn eqnFitter = " + equationFitter);
//            System.out.println("fitEqn resNums = " + resNums);
//            System.out.println("fitEqn eqnName = " + equationName);
                double[] sliderGuesses = null;
                if (sliderGuessCheckBox.isSelected()) {
                    sliderGuesses = simControls.sliderGuess(equationName, map);
                }
                CoMDOptions options = new CoMDOptions(true);
                fitResult = equationFitter.doFit(equationName, sliderGuesses, options);
                updateAfterFit(fitResult);
            }
        } catch (NullPointerException npE2) {
            npE2.printStackTrace();
            Alert alert = new Alert(Alert.AlertType.ERROR);
            alert.setContentText("Error: Residue must be selected when fitting");
            alert.showAndWait();
            return;
        }
    }

    public void updateAfterFit(FitResult fitResult) {
        List<GUIPlotEquation> equations = new ArrayList<>();
        int nCurves = fitResult.getNCurves();
        for (int iCurve = 0; iCurve < nCurves; iCurve++) {
            CurveFit curveFit = fitResult.getCurveFit(iCurve);
            List<ParValueInterface> parValues = curveFit.getParValues();
            equationChoice.getSelectionModel().select(fitResult.getEquationName());
            String aic = String.format("%.2f", fitResult.getAicc());
            String rms = String.format("%.3f", fitResult.getRms());
            String rChiSq = String.format("%.2f", fitResult.getRChiSq());
            aicLabel.setText(aic);
            rmsLabel.setText(rms);
            rChiSqLabel.setText(rChiSq);
            updateTableWithPars(parValues);
            simControls.updateSliders(parValues, fitResult.getEquationName());
            String equationName = fitResult.getEquationName(); //equationSelector.getValue();
            String expType = getFittingMode();
            //System.out.println("Fit button residueProperties = " + residueProperties);
            //System.out.println("Fit button expData = " + residueProps.getExperimentData("cest"));
            Optional<Experiment> optionalData = Optional.empty();
            if (hasExperimentSet()) {
                optionalData = getCurrentExperimentSet().getExperimentData().stream().findFirst();
            }

            if (optionalData.isPresent() && optionalData.get().getExtras().size() > 0) {
                for (Experiment expData : getCurrentExperimentSet().getExperimentData()) {
                    double[] pars = curveFit.getEquation().getPars(); //pars = getPars(equationName);
                    double[] errs = curveFit.getEquation().getErrs(); //double[] errs = new double[pars.length];
                    double[] extras = new double[3];
                    extras[0] = expData.getNucleusField();
                    extras[1] = expData.getExtras().get(0);
                    extras[2] = expData.getExtras().get(1);
//                    System.out.println("Fit button expData extras size = " + expData.getExtras().size() + " extra[1] = " + extras[1]);
                    GUIPlotEquation plotEquation = new GUIPlotEquation(expType, equationName, pars, errs, extras);
//                    for (int i = 0; i < extras.length; i++) {
//                        System.out.println(iCurve + " " + i + " extra " + extras[i]);
//                    }
//                    for (int i = 0; i < pars.length; i++) {
//                        System.out.println(iCurve + " " + i + " pars " + pars[i]);
//                    }
                    //equationCopy.setExtra(extras);

                    equations.add(plotEquation);
                }
//                ExperimentData expData = optionalData.get();
//                double[] pars = curveFit.getEquation().getPars(); //pars = getPars(equationName);
//                double[] errs = curveFit.getEquation().getErrs(); //double[] errs = new double[pars.length];
//                double[] extras = new double[3];
//                for (int j = 0; j < expData.getExtras().size() / 2; j++) {
//                    extras[0] = expData.getB0Field();
//                    extras[1] = expData.getExtras().get(2 * j);
//                    extras[2] = expData.getExtras().get(2 * j + 1);
////                    System.out.println("Fit button expData extras size = " + expData.getExtras().size() + " extra[1] = " + extras[1]);
//                    PlotEquation plotEquation = new PlotEquation(equationName, pars, errs, extras);
//                    for (int i = 0; i < extras.length; i++) {
//                        System.out.println(iCurve + " " + i + " extra " + extras[i]);
//                    }
//                    for (int i = 0; i < pars.length; i++) {
//                        System.out.println(iCurve + " " + i + " pars " + pars[i]);
//                    }
//                    //equationCopy.setExtra(extras);
//
//                    equations.add(plotEquation);
//                }
            } else {
                double[] pars = curveFit.getEquation().getPars(); //pars = getPars(equationName);
                double[] errs = curveFit.getEquation().getErrs(); //double[] errs = new double[pars.length];
                double[] extras = curveFit.getEquation().getExtras();
                double[] simExtras = simControls.getExtras();
//                extras[0] = CPMGFit.REF_FIELD;
//System.out.println("extras " + extras[0]);
                if (simExtras.length > 1) {
                    extras = new double[simExtras.length + 1];
                    extras[0] = CoMDPreferences.getRefField() * simControls.getNucleus().getFreqRatio();

                    for (int i = 0; i < simExtras.length; i++) {
                        //                    System.out.println("simextras " + i + " " + simExtras[i]);
                        extras[i + 1] = simExtras[i];
                    }
                }
                //for (int i = 0; i < extras.length; i++) {
                //System.out.println(iCurve + " " + i + " extra " + extras[i]);
                //}
                //for (int i = 0; i < simExtras.length; i++) {
                //System.out.println(iCurve + " " + i + " simExtras " + simExtras[i]);
                //}
                //for (int i = 0; i < pars.length; i++) {
                //System.out.println(iCurve + " " + i + " pars " + pars[i]);
                //}
                GUIPlotEquation plotEquation = new GUIPlotEquation(expType, equationName, pars, errs, extras);

                //equationCopy.setExtra(extras);
                //System.out.println("Fit button extras size = " + extras.length + " extra[0] = " + extras[0]);
                equations.add(plotEquation);

            }
//            double[] pars = curveFit.getEquation().getPars();
//            double[] errs = curveFit.getEquation().getErrs();
//            //System.out.println("updateAfterFit got called with curvefit.getEquation().getExtras() length = "+curveFit.getEquation().getExtras().length);
//            PlotEquation plotEquation = new PlotEquation(fitResult.getEquationName(), pars, errs, curveFit.getEquation().getExtras());
//            equations.add(plotEquation);
        }
        showEquations(equations);
    }

    @FXML
    public void fitResidues(ActionEvent event) {
//        if (getFittingMode().equals("cest")) {
//            ChooseCESTFitEquations.allRes = true;
//            ChooseCESTFitEquations.create();
//        } else {
        fitResult = null;
        if (hasExperimentSet()) {
            residueFitter.fitResidues(getCurrentExperimentSet());
        }
//        }
    }

    @FXML
    public void fitGroupResidues(ActionEvent event) {
//        if (getFittingMode().equals("cest")) {
//            ChooseCESTFitEquations.allRes = false;
//            ChooseCESTFitEquations.create();
//        } else {
        if (hasExperimentSet()) {

            fitResult = null;
            List<List<String>> allResidues = new ArrayList<>();
            List<String> groupResidues = new ArrayList<>();
            fittingResidues.clear();
            fittingResidues.addAll(ResidueChart.selectedResidues);
            groupResidues.addAll(ResidueChart.selectedResidues);
            if (!groupResidues.isEmpty()) {
                allResidues.add(groupResidues);
                residueFitter.fitResidues(getCurrentExperimentSet(), allResidues);
            }
        }
//        }
    }

    public void refreshFit() {
        ResidueChart chart = getActiveChart();
        ResidueChart.selectedResidues.clear();
        ResidueChart.selectedResidues.addAll(fittingResidues);
        refreshResidueCharts();
        chart.showInfo();
        makeAxisMenu();
    }

    @FXML
    public void haltFit(ActionEvent event) {
        residueFitter.haltFit();
        makeAxisMenu();
    }

    @FXML
    public void saveParameters(ActionEvent event) {
        FileChooser fileChooser = new FileChooser();
        fileChooser.setTitle("Save Parameter File");
        File file = fileChooser.showSaveDialog(MainApp.primaryStage);
        if (file != null) {
            DataIO.saveResultsFile(file.getAbsolutePath(), getCurrentExperimentSet(), false);
        }
    }

    @FXML
    public void saveParametersSTAR(ActionEvent event) throws IOException, InvalidMoleculeException, ParseException, InvalidPeakException {
        FileChooser fileChooser = new FileChooser();
        fileChooser.setTitle("Save STAR File");
        File file = fileChooser.showSaveDialog(MainApp.primaryStage);
        if (file != null) {
            DataIO.writeSTAR3File(file.getAbsolutePath());
            System.out.println("wrote " + file.getAbsolutePath());
        }
    }

    @FXML
    public void addT1Results(ActionEvent event) {
        MoleculeBase mol = MoleculeFactory.getActive();
        String alertText = "Add T1 results to map?";
        if (mol != null) {
            alertText = "Add T1 results to molecule " + mol.getName() + "?";
        }
        Alert alert = new Alert(Alert.AlertType.CONFIRMATION, alertText);
        Optional<ButtonType> response = alert.showAndWait();
        if (response.isPresent() && response.get().getText().equals("OK")) {
            DataIO.addRelaxationFitResults(getCurrentExperimentSet(), relaxTypes.T1);
        }
    }

    @FXML
    public void addT2Results(ActionEvent event) {
        MoleculeBase mol = MoleculeFactory.getActive();
        String alertText = "Add T2 results to map?";
        if (mol != null) {
            alertText = "Add T2 results to molecule " + mol.getName() + "?";
        }
        Alert alert = new Alert(Alert.AlertType.CONFIRMATION, alertText);
        Optional<ButtonType> response = alert.showAndWait();
        if (response.isPresent() && response.get().getText().equals("OK")) {
            DataIO.addRelaxationFitResults(getCurrentExperimentSet(), relaxTypes.T2);
        }
    }

    public Double updateFitProgress(Double f) {
        if (hasExperimentSet()) {
            ResidueChart chart = getActiveChart();
            String seriesName = chart.currentSeriesName;
            String[] sParts;
            if (seriesName.length() == 0) {
                sParts = new String[4];
                sParts[0] = getCurrentExperimentSet().getName();
                sParts[1] = "best";
                sParts[2] = "0:0:0";
                sParts[3] = "RMS";
            } else {
                sParts = seriesName.split("\\|");
            }
            if (Platform.isFxApplicationThread()) {
                clearChart();
                statusBar.setProgress(f);
                setYAxisType(getCurrentExperimentSet().getExpMode(), sParts[0], sParts[1], sParts[2], sParts[3]);
                setCurrentExperimentSet(ChartUtil.getResidueProperty(getCurrentExperimentSet().getName()));

            } else {
                Platform.runLater(() -> {
                    clearChart();
                    setYAxisType(getCurrentExperimentSet().getExpMode(), sParts[0], sParts[1], sParts[2], sParts[3]);
                    statusBar.setProgress(f);
                    setCurrentExperimentSet(ChartUtil.getResidueProperty(getCurrentExperimentSet().getName()));

                });
            }
        }
        return null;

    }

    public void processingDone() {
        makeAxisMenu();
    }

    public Double updateStatus(ProcessingStatus status) {
        if (Platform.isFxApplicationThread()) {
            updateStatusNow(status);
        } else {
            Platform.runLater(() -> {
                updateStatusNow(status);
            });

        }
        return null;
    }

    public Double updateStatusNow(ProcessingStatus status) {
        String s = status.getStatus();
        if (s == null) {
            statusBar.setText("");
        } else {
            statusBar.setText(s);
            if (s.equals("Done")) {
                refreshFit();
            }
        }
        if (status.isOk()) {
            statusCircle.setFill(Color.GREEN);
        } else {
            if (status.getThrowable() != null) {
                status.getThrowable().printStackTrace();
            }
            statusCircle.setFill(Color.RED);
        }
        statusBar.setProgress(0.0);
        return null;
    }

    public void saveBarChart() throws IOException {
        FileChooser chooser = new FileChooser();
        chooser.setInitialFileName("barchart.png");
        File file = chooser.showSaveDialog(MainApp.primaryStage);
        if (file != null) {
            snapit(chartBox, file);
        }
    }

    public void saveXYChart() throws IOException {
        FileChooser chooser = new FileChooser();
        chooser.setInitialFileName("cpmgchart.png");
        File file = chooser.showSaveDialog(MainApp.primaryStage);
        if (file != null) {
            snapit(chartPane, file);
        }
    }

    public void exportExecutable(String fileSuggestion, String exportType, boolean saveBar) throws ScriptException {
        FileChooser chooser = new FileChooser();
        chooser.setInitialFileName(fileSuggestion);
        File file = chooser.showSaveDialog(MainApp.primaryStage);
        if (file != null) {
            String filePath = file.getAbsolutePath();

            Map<String, Object>[] plottedData = xychart.getPlottedData();
            Map<String, Object> graphData = null;
            graphData = xychart.getGraphData();
            graphData.put("file", filePath);
            graphData.put("exportType", exportType);
            ArrayList<Object> barChartData;
            if (!"grace".equals(exportType) && saveBar) {
                ObservableList<Node> barNodes = chartBox.getChildren();
                barChartData = new ArrayList<>(barNodes.size());
                barNodes.forEach(node -> {
//                    if (node instanceof XYBarChart) {
//                        XYBarChart barChart = (XYBarChart) node;
//                        String barChartName = barChart.toString();
//                        HashMap<String, Object> chartData = barChart.getChartData();
//                        barChartData.add(chartData);
//                    }
                });
            } else {
                barChartData = new ArrayList<>(0);
            }

            MainApp.engine.put("data", plottedData);
            MainApp.engine.put("configData", graphData);
            MainApp.engine.put("barChartData", barChartData);
            MainApp.engine.eval("writer = Writer(configData,data,barChartData)");
            MainApp.engine.eval("writer.writeFromExportData()");
        }
    }

    public void saveGraceFile() throws IOException, ScriptException {
        exportExecutable("ASCII.agr", "grace", false);
    }

    public void savePythonFile() throws IOException, ScriptException {
        exportExecutable("graph.py", "python", false);
    }

    public void saveRFile() throws IOException, ScriptException {
        exportExecutable("graph.r", "r", false);
    }

    public void saveBarToPythonFile() throws IOException, ScriptException {
        exportExecutable("graph.py", "python", true);
    }

    public void saveBarToRFile() throws IOException, ScriptException {
        exportExecutable("graph.r", "r", true);
    }

    public void snapit(Node node, File file) throws IOException {
        double scale = 4.0;
        final Bounds bounds = node.getLayoutBounds();
        final WritableImage image = new WritableImage(
                (int) Math.round(bounds.getWidth() * scale),
                (int) Math.round(bounds.getHeight() * scale));
        final SnapshotParameters spa = new SnapshotParameters();
        spa.setTransform(javafx.scene.transform.Transform.scale(scale, scale));
        node.snapshot(spa, image);
        ImageIO.write(SwingFXUtils.fromFXImage(image, null), "png", file);
    }

    public FitResult getFitResult() {
        return fitResult;
    }

    public ExperimentResult getResidueInfo() {
        return chartInfo.getResult();
    }

    public String getEquationName() {
        return chartInfo.getEquationName();
    }

    public String getFittingMode() {
        String fitMode = "cpmg";
        if (!hasExperimentSet()) {
            String simMode = getSimMode();
            if (simMode != null) {
                fitMode = simMode;
            }
        } else {
            fitMode = getCurrentExperimentSet().getExpMode();
            if (fitMode.equalsIgnoreCase("t1") || fitMode.equalsIgnoreCase("t2")) {
                fitMode = "exp";
            }
        }
        return fitMode.toLowerCase();
    }

    public EquationFitter getFitter() {
        CoMDOptions options = new CoMDOptions(true);
        if (getFittingMode().equals("exp")) {
            return new ExpFitter(options);
        } else if (getFittingMode().equals("cpmg")) {
            return new CPMGFitter(options);
        } else if (getFittingMode().equals("cest")) {
            return new CESTFitter(options);
        } else if (getFittingMode().equals("r1rho")) {
            return new R1RhoFitter(options);
        }
        return null;
    }

    public String[] getParNames(String equationName) {
        String fitMode = getFittingMode();
        EquationType type = ResidueFitter.getEquationType(fitMode, equationName);
        return type.getParNames();
    }

    public String[] getParTypes(String mode) {
        String[] cpmgTypes = {"R2", "dPPMmin", "Kex", "pA", "dPPM", "RMS", "AIC", "Equation"};
        String[] expTypes = {"A", "R", "C", "RMS", "AIC", "Equation"};
        String[] t1Types = {"R1"};
        String[] t2Types = {"R2"};
        String[] sTypes = {"S2", "Rex"};
        String[] cestTypes = {"kex", "pb", "deltaA0", "deltaB0", "R1A", "R1B", "R2A", "R2B", "RMS", "AIC", "Equation"};
        String[] r1rhoTypes = {"kex", "pb", "deltaA0", "deltaB0", "R1A", "R1B", "R2A", "R2B", "RMS", "AIC", "Equation"};
        String[] nullTypes = {"RMS", "AIC", "Equation"};
        String[] noeTypes = {"NOE"};
        if (mode.equals("exp")) {
            return expTypes;
        } else if (mode.equals("cpmg")) {
            return cpmgTypes;
        } else if (mode.equals("cest")) {
            return cestTypes;
        } else if (mode.equals("r1rho")) {
            return r1rhoTypes;
        } else if (mode.equals("t1")) {
            return expTypes;
        } else if (mode.equals("t2")) {
            return expTypes;
        } else if (mode.equals("s2")) {
            return sTypes;
        } else if (mode.equals("noe")) {
            return noeTypes;
        }
        return nullTypes;
    }

    @FXML
    void setBestEquation(ActionEvent e) {
        for (String resNum : chartInfo.currentResidues) {
            ExperimentResult resInfo = ChartUtil.getResInfo(chartInfo.mapName, String.valueOf(resNum));
            if (resInfo != null) {
                String equationName = equationChoice.getValue();
                resInfo.setBestEquationName(equationName);
            }
        }
        updateFitProgress(1.0);

    }

    void equationAction() {
        if (equationChoice.getUserData() == null) {
            String equationName = equationChoice.getValue();
            if (!chartInfo.currentStates.isEmpty() && equationName != null) {
                // copy it so it doesn't get cleared by clear call in updateTableWithPars
                updateTableWithPars(chartInfo, false);
                showInfo(equationName);
            }
        }
    }

    public String getParametersEquation() {
        return equationChoice.getValue();
    }

    public String getSimMode() {
        String simSelected = simChoice.getSelectionModel().getSelectedItem();
        if (simSelected == null) {
            return "cpmg";
        }
        String simMode = simSelected.substring(9, simSelected.length()).toLowerCase();
        return simMode;
    }

    void simAction() {
        clearProject(false);
        getSimMode();
        setSimControls();
        updateXYChartLabels();
        genDataSDevTextField.setText("");
        simControls.simSliderAction("");
    }

    @FXML
    void clearProject(ActionEvent event) {
        clearProject(true);
    }

    void clearProject(boolean clearXY) {
        chartInfo.clear();
        setCurrentExperimentSet(null);
        fitResult = null;
        ChartUtil.clearResidueProperties();
        if (clearXY) {
            xychart.clear();
        } else {
            xychart.getData().clear();
        }
        clearSecondaryStructure();
        barCharts.remove(activeChart);
        chartBox.getChildren().remove(0, chartBox.getChildren().size());
        chartBox.getChildren().add(barPlotCanvas);
        addChart();
    }

    void showInfo(String equationName) {
        if (hasExperimentSet()) {
            showInfo(chartInfo, xychart);
        }
    }

    void showInfo(ChartInfo chartInfo, PlotData plotData) {
        ArrayList<GUIPlotEquation> equations = new ArrayList<>();
        ObservableList<DataSeries> allData = FXCollections.observableArrayList();
        List<ExperimentData> experimentalDataSets = new ArrayList<>();
        List<int[]> allStates = new ArrayList<>();
        boolean calcScale = scalePlot.isSelected();
        System.out.println("showinfo " + chartInfo.hasExperiments() + " " + chartInfo.hasResidues());
        if (chartInfo.hasExperiments() && chartInfo.hasResidues()) {
            //double maxY = getMaxY(experimentSet, equationName, mapName, state, residues) / 100.0;
            //System.out.println("max Y " + maxY);
            int iSeries = 0;
            for (Experiment expData : chartInfo.getExperiments().getExperimentData()) {
                if (!ExperimentSet.matchStateString(chartInfo.state, expData.getState())) {
                    continue;
                }
                String expName = expData.getName();
                for (String resNum : chartInfo.getResidues()) {
                    if (expData.getResidueData(resNum) != null) {
//                        System.out.println(expData.getResidueData(resNum));
                        experimentalDataSets.add(expData.getResidueData(resNum));
                        DataSeries series = ChartUtil.getMapData(chartInfo.mapName, expName, resNum);
                        series.setStroke(PlotData.colors[iSeries % 8]);
                        series.setFill(PlotData.colors[iSeries % 8]);
                        allData.add(series);
                        GUIPlotEquation equation = ChartUtil.getEquation(expData,
                                chartInfo.mapName, resNum, chartInfo.equationName, expData.getState(),
                                expData.getNucleusField());
                        double maxY = 1.0;
                        if (equation != null) {
                            equation.setColor(PlotData.colors[iSeries % 8]);
                            if (calcScale) {
                                maxY = equation.calculate(equation.getMinX()) / 100.0;
                            }
                            equation.setScaleValue(maxY);
                            equations.add(equation);
                        } else if (calcScale) {
                            maxY = series.getValues().stream().mapToDouble(XYValue::getYValue).max().getAsDouble() / 100.0;
                        }
                        series.setScale(maxY);
                        iSeries++;
                    }
                }

                int[] states = chartInfo.currentExperimentSet.getStateIndices(0, expData);
                allStates.add(states);
            }
        }
        chartInfo.setStates(allStates);
        updateTable(experimentalDataSets);
        if (chartInfo.getResidues() != null) {
            updateTableWithPars(chartInfo, true);
            updateEquation(chartInfo.mapName, chartInfo.getResidues(), chartInfo.equationName);
        }
        plotData.setData(allData);
        setBounds();
        plotData.autoScale(false);
        plotData.setEquations(equations);
    }

    public double getMaxY(ExperimentSet experimentSet, String equationName, String mapName, String state, String[] residues) {
        double maxValue = Double.NEGATIVE_INFINITY;
        if ((experimentSet != null) && (residues != null)) {
            for (Experiment expData : experimentSet.getExperimentData()) {
                if (!ExperimentSet.matchStateString(state, expData.getState())) {
                    continue;
                }
                String expName = expData.getName();
                for (String resNum : residues) {

                    GUIPlotEquation equation = ChartUtil.getEquation(expData,
                            mapName, resNum, equationName, expData.getState(),
                            expData.getNucleusField());
                    if (equation != null) {
                        double minX = equation.getMinX();
                        double valueY = equation.calculate(minX);
                        if (valueY > maxValue) {
                            maxValue = valueY;
                        }
                    } else {
                        DataSeries series = ChartUtil.getMapData(mapName, expName, resNum);
                        double valueY = series.getValues().stream().mapToDouble(XYValue::getYValue).max().getAsDouble();
                        if (valueY > maxValue) {
                            maxValue = valueY;
                        }

                    }
                }
            }
        }
        return maxValue;
    }

    public void setSimData(EquationFitter equationFitter) {
        ArrayList<Double> xValues = getSimXData();
        ArrayList<Double> yValues = getSimYData();
        ArrayList<Double> errValues = getSimErrData();
        double[] extras = simControls.getExtras();
        double[] fieldVals = new double[yValues.size()];
        ArrayList[] allXValues = new ArrayList[extras.length + 1];
        allXValues[0] = xValues;
        ArrayList<Double> fieldValues = new ArrayList<>();
        for (int i = 0; i < yValues.size(); i++) {
            fieldVals[i] = CoMDPreferences.getRefField() * simControls.getNucleus().getFreqRatio();
            fieldValues.add(fieldVals[i]);
        }
        System.out.println("xval len " + xValues.size());
        for (int j = 0; j < extras.length; j++) {
            System.out.println("set sim " + j + " " + extras[j]);
            ArrayList<Double> xValuesEx = new ArrayList<>();
            for (int i = 0; i < yValues.size(); i++) {
                xValuesEx.add(extras[j]);
            }
            allXValues[1 + j] = xValuesEx;
        }
        equationFitter.getFitModel().setFieldValues(fieldVals);
        equationFitter.setData(allXValues, yValues, errValues, fieldValues);
    }

    public void guessSimData() {
        EquationFitter equationFitter = getFitter();
        setSimData(equationFitter);
        String equationName = simControls.getEquation();
        List<ParValueInterface> guesses = equationFitter.guessPars(equationName);
        if (guesses != null) {
            simControls.updateSliders(guesses, equationName);
            simControls.simSliderAction("");
        }
    }

    public void fitSimData() {
        EquationFitter equationFitter = getFitter();
        setSimData(equationFitter);
        String equationName = simControls.getEquation();
        equationFitter.setupFit(equationName);
        int[][] map = equationFitter.getFitModel().getMap();
        double[] sliderGuesses = null;
        if (sliderGuessCheckBox.isSelected()) {
            sliderGuesses = simControls.sliderGuess(equationName, map);
        }
        CoMDOptions options = new CoMDOptions(true);
        fitResult = equationFitter.doFit(equationName, sliderGuesses, options);
        updateAfterFit(fitResult);
    }

    ArrayList<Double> getSimXData() {
        ObservableList<DataSeries> allData = xychart.getData();
        ArrayList<Double> xValues = new ArrayList<>();
        for (DataSeries series : allData) {
            for (XYValue dataPoint : series.getData()) {
                double x = dataPoint.getXValue();
                xValues.add(x);
            }
        }
        return xValues;
    }

    ArrayList<Double> getSimYData() {
        ObservableList<DataSeries> allData = xychart.getData();
        ArrayList<Double> yValues = new ArrayList<>();
        for (DataSeries series : allData) {
            for (XYValue dataPoint : series.getData()) {
                double y = dataPoint.getYValue();
                yValues.add(y);
            }
        }
        return yValues;
    }

    ArrayList<Double> getSimErrData() {
        ObservableList<DataSeries> allData = xychart.getData();
        ArrayList<Double> errValues = new ArrayList<>();
        for (DataSeries series : allData) {
            for (XYValue dataPoint : series.getData()) {
                double errValue = 0.0;
                if (dataPoint instanceof XYEValue) {
                    XYEValue xyeValue = (XYEValue) dataPoint;
                    errValue = xyeValue.getError();
                }
                errValues.add(errValue);
            }
        }
        return errValues;
    }

    @FXML
    void showSimData(ActionEvent e) {
        ObservableList<DataSeries> allData = FXCollections.observableArrayList();
        String equationName = simControls.getEquation();
        String simMode = getSimMode();
        EquationType eType = ResidueFitter.getEquationType(simMode, equationName);
        int[][] map = eType.makeMap(1);
        EquationFitter equationFitter = getFitter();
        equationFitter.getFitModel().setMap(map);
        double[] sliderGuesses = simControls.sliderGuess(equationName, map);
        double[] xBounds = xychart.getXBounds();
        double[] yBounds = xychart.getYBounds();
        double sdev = Math.abs(yBounds[1] - yBounds[0]) * 0.02;
        if (genDataSDevTextField.getText().equals("")) {
            genDataSDevTextField.setText(String.valueOf(sdev));
        }
        double[] xValues = equationFitter.getSimXDefaults();
        if (simMode.equals("cest") || simMode.equals("r1rho")) {
            int nPts = Integer.parseInt(genDataNPtsTextField.getText());
            double xLB = Double.parseDouble(genDataXLBTextField.getText());
            double xUB = Double.parseDouble(genDataXUBTextField.getText());
            xValues = equationFitter.getSimX(nPts, xLB, xUB);
        } else {
            String[] xvals = genDataXValTextField.getText().split(" ");
            double[] xVals = new double[xvals.length];
            for (int i = 0; i < xvals.length; i++) {
                xVals[i] = Double.parseDouble(xvals[i]);
            }
            xValues = xVals;
        }
        double fieldRef = 1.0;
        int iLine = 0;

        List<DataSeries> data = new ArrayList<>();
        DataSeries series = new DataSeries();
        series.setName("sim" + ":" + "0");
        data.add(series);
        for (PlotEquation eqn : xychart.plotEquations) {
            if (iLine == 0) {
                fieldRef = eqn.getExtra(0);
            }
//            plotEquation.equation.setFieldRef(fieldRef);
            double[] extras = eqn.getExtras();
            double[] ax = new double[extras.length];
//            System.out.println("extras " + extras.length);
            for (int j = 1; j < extras.length; j++) {
                ax[j] = extras[j];
            }
            for (int i = 0; i < xValues.length; i++) {
                double xValue = xValues[i];
                ax[0] = xValue;
                double yValue = eqn.calculate(sliderGuesses, ax, fieldRef);
                yValue += Double.parseDouble(genDataSDevTextField.getText()) * rand.nextGaussian(); //sdev * rand.nextGaussian();
//                XYValue dataPoint = new XYValue(xValue, yValue);
                XYValue dataPoint = new XYEValue(xValue, yValue, Double.parseDouble(genDataSDevTextField.getText()));
                dataPoint.setExtraValue(new Double(sdev));
                series.getData().add(dataPoint);
            }
        }
        allData.addAll(data);
        xychart.setData(allData);
//    
    }

    @FXML
    public void loadSimData(ActionEvent e) {
        ObservableList<DataSeries> allData = FXCollections.observableArrayList();
        FileChooser fileChooser = new FileChooser();
        File file = fileChooser.showOpenDialog(null);
        if (file != null) {
            try {
                List<Double>[] dataValues = DataIO.loadSimDataFile(file);
                List<DataSeries> data = new ArrayList<>();
                DataSeries series = new DataSeries();
                series.setName("sim" + ":" + "0");
                data.add(series);
                for (int i = 0; i < dataValues[0].size(); i++) {
                    XYValue dataPoint;
                    if (!dataValues[2].isEmpty()) {
                        dataPoint = new XYEValue(dataValues[0].get(i),
                                dataValues[1].get(i), dataValues[2].get(i));
                    } else {
                        dataPoint = new XYValue(dataValues[0].get(i), dataValues[1].get(i));

                    }
                    series.getData().add(dataPoint);
                }
                allData.addAll(data);
                xychart.setData(allData);
            } catch (IOException ioE) {
                ExceptionDialog dialog = new ExceptionDialog(ioE);
                dialog.showAndWait();
            }
        }
    }

    @FXML
    public void showMCplot(ActionEvent event) {
        if (bootstrapSamplePlots == null) {
            bootstrapSamplePlots = new BootstrapSamplePlots(this);
        }
        bootstrapSamplePlots.showMCplot();
    }

    @FXML
    private void showPreferences(ActionEvent event) {
        if (preferencesController == null) {
            preferencesController = PreferencesController.create(primaryStage);
        }
        if (preferencesController != null) {
            preferencesController.getStage().show();
        } else {
            System.out.println("Coudn't make controller");
        }
    }

    @FXML
    private void showConsole(ActionEvent event) {
        if (console == null) {
            console = ConsoleRedirect.create();
        }
        if (console != null) {
            console.show();
        } else {
            System.out.println("Coudn't make console");
        }
    }

    @FXML
    private void exitProgram() {
        MainApp.checkExit();
    }

    @FXML
    private void printXYChart() {
        printChart(chartPane);
    }

    @FXML
    private void printBarChart() {
        printChart(chartBox);
    }

    private void printChart(Node node) {
        PrinterJob job = PrinterJob.createPrinterJob();
        if (job != null) {
            boolean selected = job.showPrintDialog(null);
            if (selected) {
                boolean pageSetup = job.showPageSetupDialog(null);
                if (pageSetup) {
                    boolean success = job.printPage(node);
                    if (success) {
                        job.endJob();
                    }
                }
            }
        }
    }

    public File getInitialDirectory() {
        if (initialDir == null) {
            String homeDirName = System.getProperty("user.home");
            initialDir = new File(homeDirName);
        }
        return initialDir;
    }

    @FXML
    void exportSVGAction(ActionEvent event) {
        FileChooser fileChooser = new FileChooser();
        fileChooser.setTitle("Export to SVG");
        fileChooser.setInitialDirectory(getInitialDirectory());
        File selectedFile = fileChooser.showSaveDialog(null);
        if (selectedFile != null) {
            SVGGraphicsContext svgGC = new SVGGraphicsContext();
            try {
                Canvas canvas = xychart.getCanvas();
                svgGC.create(true, canvas.getWidth(), canvas.getHeight(), selectedFile.toString());
                xychart.exportVectorGraphics(svgGC);
                svgGC.saveFile();
            } catch (GraphicsIOException ex) {
                ExceptionDialog eDialog = new ExceptionDialog(ex);
                eDialog.showAndWait();
            }
        }
    }

    @FXML
    void exportBarPlotSVGAction(ActionEvent event) {
        FileChooser fileChooser = new FileChooser();
        fileChooser.setTitle("Export to SVG");
        fileChooser.setInitialDirectory(getInitialDirectory());
        File selectedFile = fileChooser.showSaveDialog(null);
        if (selectedFile != null) {
            SVGGraphicsContext svgGC = new SVGGraphicsContext();
            try {
                svgGC.create(true, barPlotCanvas.getWidth(), barPlotCanvas.getHeight(), selectedFile.toString());
                exportBarPlotSVGAction(svgGC);
                svgGC.saveFile();
            } catch (GraphicsIOException ex) {
                ExceptionDialog eDialog = new ExceptionDialog(ex);
                eDialog.showAndWait();
            }
        }
    }

    protected void exportBarPlotSVGAction(SVGGraphicsContext svgGC) throws GraphicsIOException {
        svgGC.beginPath();
        for (ResidueChart resChart : barCharts) {
            resChart.drawChart(svgGC);
        }
        if (ssPainter != null) {
            double ssHeight = ssPainter.getHeight();
            Axis axis = activeChart.xAxis;
            ssPainter.paintSS(svgGC, axis.getXOrigin(), barPlotCanvas.getHeight() - ssHeight,
                    axis.getWidth(),
                    axis.getLowerBound(), axis.getUpperBound());
        }
    }

    public boolean hasExperimentSet() {
        return (chartInfo.currentExperimentSet != null) && (chartInfo.currentExperimentSet instanceof ExperimentSet);
    }

    /**
     * @return the currentExperimentSet
     */
    public ExperimentSet getCurrentExperimentSet() {
        return chartInfo.currentExperimentSet;
    }

    /**
     * @param currentResProps the currentExperimentSet to set
     */
    public void setCurrentExperimentSet(ExperimentSet currentResProps) {
        this.chartInfo.currentExperimentSet = currentResProps;
    }

}

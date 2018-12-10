package org.comdnmr.fit.gui;

import de.jensd.fx.glyphs.GlyphsDude;
import de.jensd.fx.glyphs.fontawesome.FontAwesomeIcon;
import java.io.File;
import java.io.IOException;
import javafx.fxml.FXML;
import javafx.fxml.Initializable;
import java.util.ResourceBundle;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import javafx.application.Platform;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.embed.swing.SwingFXUtils;
import javafx.event.ActionEvent;
import javafx.event.Event;
import javafx.geometry.Bounds;
import javafx.geometry.Side;
import javafx.scene.Node;
import javafx.scene.SnapshotParameters;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.XYChart;
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
import javafx.scene.layout.Priority;
import javafx.scene.layout.VBox;
import javafx.scene.paint.Color;
import javafx.scene.shape.Circle;
import javafx.stage.FileChooser;
import javafx.stage.Stage;
import javax.imageio.ImageIO;
import javax.script.ScriptException;
import org.comdnmr.fit.calc.CPMGFit;
import org.comdnmr.fit.calc.CPMGFitResult;
import org.comdnmr.fit.calc.CurveFit;
import org.comdnmr.fit.calc.DataIO;
import org.comdnmr.fit.calc.EquationFitter;
import org.comdnmr.fit.calc.EquationType;
import org.comdnmr.fit.calc.ExpFit;
import org.comdnmr.fit.calc.ExperimentData;
import org.comdnmr.fit.calc.ResidueData;
import org.controlsfx.control.PropertySheet;
import org.comdnmr.fit.calc.ParValueInterface;
import org.comdnmr.fit.calc.PlotEquation;
import org.comdnmr.fit.calc.ProcessingStatus;
import org.comdnmr.fit.calc.ResidueFitter;
import org.comdnmr.fit.calc.ResidueInfo;
import org.comdnmr.fit.calc.ResidueProperties;
import org.controlsfx.control.StatusBar;
import org.comdnmr.fit.calc.CESTFit;
import java.text.DecimalFormat;
import java.util.Collections;
import java.util.HashMap;
import java.util.Optional;
import javafx.beans.property.SimpleStringProperty;
import java.util.Random;
import javafx.scene.control.Alert;
import javafx.scene.control.SplitPane;
import javafx.scene.control.Label;
import javafx.scene.control.TextField;
import javafx.scene.control.Button;
import javafx.scene.control.ContentDisplay;
import javafx.scene.control.ToolBar;
import org.comdnmr.fit.calc.FitModel;
import static org.comdnmr.fit.gui.MainApp.preferencesController;
import static org.comdnmr.fit.gui.MainApp.console;
import static org.comdnmr.fit.gui.MainApp.primaryStage;
import org.controlsfx.dialog.ExceptionDialog;

public class PyController implements Initializable {

    public static PyController mainController;
    XYBarChart activeChart;
    @FXML
    XYBarChart xyBarChart0;

    SSRegion activeSSregion;
    @FXML
    SSRegion ssregion;

    @FXML
    PlotData xychart;

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
    VBox chartBox;
    @FXML
    StatusBar statusBar;
    Circle statusCircle;
    @FXML
    CheckBox absValueModeCheckBox;
    @FXML
    CheckBox nonParBootStrapCheckBox;
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
    Button setBoundsButton;

    EquationControls simControls;

    @FXML
    SplitPane splitPane;

    BootstrapSamplePlots bootstrapSamplePlots = null;
    InputDataInterface inputDataInterface = null;

    //final ContextMenu axisMenu = new ContextMenu();
    static double defaultField = 500.0;

    ResidueInfo currentResInfo = null;
    ResidueProperties currentResProps = null;
    ResidueFitter residueFitter;
    String currentMapName = "";
    String[] currentResidues;
    String currentState = "";
    String currentEquationName;
    List<int[]> currentStates = new ArrayList<>();

    boolean simulate = true;

    CPMGFitResult fitResult;

    @FXML
    ToolBar navigatorToolBar = new ToolBar();

    static Random rand = new Random();

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

    @Override
    public void initialize(URL url, ResourceBundle rb) {
        mainController = this;
        activeChart = xyBarChart0;
        activeSSregion = ssregion;
        //propertySheet.setPropertyEditorFactory(new NvFxPropertyEditorFactory());
//        propertySheet.setMode(PropertySheet.Mode.NAME);
//        propertySheet.setModeSwitcherVisible(false);
//        propertySheet.setSearchBoxVisible(false);
        makeAxisMenu();
//        simControls = new CPMGControls();
        if (getFittingMode().equals("cpmg")) {
            simControls = new CPMGControls();
            equationChoice.getItems().clear();
            equationChoice.getItems().add("+");
            equationChoice.getItems().addAll(CPMGFit.getEquationNames());
            equationChoice.setValue(CPMGFit.getEquationNames().get(0));
            xLowerBoundTextField.setText("0.0");
            xUpperBoundTextField.setText("1100.0");
            yLowerBoundTextField.setText("0.0");
            yUpperBoundTextField.setText("65.0");
            xTickTextField.setText("100.0");
            yTickTextField.setText("5.0");
        } else if (getFittingMode().equals("cest")) {
            simControls = new CESTControls();
            equationChoice.getItems().clear();
            equationChoice.getItems().add("+");
            equationChoice.getItems().addAll(CESTFit.getEquationNames());
            equationChoice.setValue(CESTFit.getEquationNames().get(0));
            xLowerBoundTextField.setText("-6000.0");
            xUpperBoundTextField.setText("6000.0");
            if (currentResProps.getExperimentData() != null) {
                double[] xVals = currentResProps.getExperimentData().stream().findFirst().get().getXVals();
                xLowerBoundTextField.setText(String.valueOf(Math.round(xVals[1]/1000)*1000));
                xUpperBoundTextField.setText(String.valueOf(Math.round(xVals[xVals.length-1]/1000)*1000));
            }
            yLowerBoundTextField.setText("0.0");
            yUpperBoundTextField.setText("1.0");
            xTickTextField.setText("2000.0");
            yTickTextField.setText("0.25");
        } else if (getFittingMode().equals("exp")) {
            simControls = new ExpControls();
            equationChoice.getItems().clear();
            equationChoice.getItems().add("+");
            equationChoice.getItems().addAll(ExpFit.getEquationNames());
            equationChoice.setValue(ExpFit.getEquationNames().get(0));
            xLowerBoundTextField.setText("0.0");
            xUpperBoundTextField.setText("1.25");
            yLowerBoundTextField.setText("0.0");
            yUpperBoundTextField.setText("200.0");
            xTickTextField.setText("0.25");
            yTickTextField.setText("50.0");
        } else {
            System.out.println("Error: no fitting mode selected.");
        }
        VBox vBox = simControls.makeControls(mainController);
        simPane.centerProperty().set(vBox);
        residueFitter = new ResidueFitter(this::updateFitProgress, this::updateStatus);
        statusCircle = new Circle(10);
        statusBar.getLeftItems().add(statusCircle);
        equationChoice.valueProperty().addListener(e -> {
            equationAction();
        });

        simChoice.getItems().add("Simulate CPMG");
        simChoice.getItems().add("Simulate EXP");
        simChoice.getItems().add("Simulate CEST");
        simChoice.setValue("Simulate CPMG");
        simChoice.valueProperty().addListener(s -> {
            simAction();
        });

        splitPane.setDividerPositions(0.4, 0.7);

        setBoundsButton.setOnAction(this::setBounds);

        initResidueNavigator();
        calcErrorsCheckBox.selectedProperty().addListener(e -> FitModel.setCalcError(calcErrorsCheckBox.isSelected()));
        calcErrorsCheckBox.setSelected(true);

//        mainController.setOnHidden(e -> Platform.exit());
    }

    public void setControls() {
        boolean update = false;
        if (getFittingMode().equals("cpmg") && !(simControls instanceof CPMGControls)) {
            simControls = new CPMGControls();
            update = true;
            equationChoice.getItems().clear();
            equationChoice.getItems().add("+");
            equationChoice.getItems().addAll(CPMGFit.getEquationNames());
            equationChoice.setValue(CPMGFit.getEquationNames().get(0));

        } else if (getFittingMode().equals("exp") && !(simControls instanceof ExpControls)) {
            simControls = new ExpControls();
            equationChoice.getItems().clear();
            equationChoice.getItems().add("+");
            equationChoice.getItems().addAll(ExpFit.getEquationNames());
            equationChoice.setValue(ExpFit.getEquationNames().get(0));
            update = true;
        } else if (getFittingMode().equals("cest") && !(simControls instanceof CESTControls)) {
            simControls = new CESTControls();
            equationChoice.getItems().clear();
            equationChoice.getItems().add("+");
            equationChoice.getItems().addAll(CESTFit.getEquationNames());
            equationChoice.setValue(CESTFit.getEquationNames().get(0));
            update = true;
        }
        if (update) {
            VBox vBox = simControls.makeControls(mainController);
            simPane.centerProperty().set(vBox);
        }
        updateXYChartLabels();
    }

    public void setSimControls() {
        boolean update = false;
        if (getSimMode().equals("cpmg") && !(simControls instanceof CPMGControls)) {
            simControls = new CPMGControls();
            update = true;
            equationChoice.getItems().clear();
            equationChoice.getItems().add("+");
            equationChoice.getItems().addAll(CPMGFit.getEquationNames());
            equationChoice.setValue(CPMGFit.getEquationNames().get(0));

        } else if (getSimMode().equals("exp") && !(simControls instanceof ExpControls)) {
            simControls = new ExpControls();
            equationChoice.getItems().clear();
            equationChoice.getItems().add("+");
            equationChoice.getItems().addAll(ExpFit.getEquationNames());
            equationChoice.setValue(ExpFit.getEquationNames().get(0));
            update = true;
        } else if (getSimMode().equals("cest") && !(simControls instanceof CESTControls)) {
            simControls = new CESTControls();
            ((CESTControls) simControls).updateDeltaLimits();
            equationChoice.getItems().clear();
            equationChoice.getItems().add("+");
            equationChoice.getItems().addAll(CESTFit.getEquationNames());
            equationChoice.setValue(CESTFit.getEquationNames().get(0));
            update = true;
        }
        if (update) {
            VBox vBox = simControls.makeControls(mainController);
            simPane.centerProperty().set(vBox);
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
        List<ResidueInfo> resInfo = currentResProps.getResidueValues();
        List resNums = new ArrayList<>();
        for (int i = 0; i < resInfo.size(); i++) {
            resNums.add(resInfo.get(i).getResNum());
        }
        Collections.sort(resNums);
        if (currentResInfo != null) {
            int resIndex = resNums.indexOf(currentResInfo.getResNum());
            resIndex--;
            if (resIndex <= 0) {
                resIndex = 0;
            }
            int res = (int) resNums.get(resIndex);
            XYBarChart chart = getActiveChart();
            chart.showInfo(chart.currentSeriesName, 0, res, false);
        }
    }

    public void firstResidue(ActionEvent event) {
        if (currentResidues != null) {
            int res = ChartUtil.minRes;
            XYBarChart chart = getActiveChart();
            chart.showInfo(chart.currentSeriesName, 0, res, false);
        }
    }

    public void nextResidue(ActionEvent event) {
        List<ResidueInfo> resInfo = currentResProps.getResidueValues();
        List resNums = new ArrayList<>();
        for (int i = 0; i < resInfo.size(); i++) {
            resNums.add(resInfo.get(i).getResNum());
        }
        Collections.sort(resNums);
        if (currentResInfo != null) {
            int resIndex = resNums.indexOf(currentResInfo.getResNum());
            resIndex++;
            if (resIndex >= resNums.size()) {
                resIndex = resNums.size() - 1;
            }
            int res = (int) resNums.get(resIndex);
            XYBarChart chart = getActiveChart();
            chart.showInfo(chart.currentSeriesName, 0, res, false);
        }
    }

    public void lastResidue(ActionEvent event) {
        if (currentResidues != null) {
            int res = ChartUtil.maxRes;
            XYBarChart chart = getActiveChart();
            chart.showInfo(chart.currentSeriesName, 0, res, false);
        }
    }

    public void loadParameterFile(Event e) {
        FileChooser fileChooser = new FileChooser();
        fileChooser.getExtensionFilters().addAll(new FileChooser.ExtensionFilter("Yaml File", "*.yaml", "*.yml"));
        Stage stage = MainApp.primaryStage;
        File file = fileChooser.showOpenDialog(stage);
        if (file != null) {
            if (activeChart != null) {
                clearChart();
                currentResidues = null;
                simulate = false;
                fitResult = null;
            }
            ChartUtil.loadParameters(file.toString());
        }
    }

    @FXML
    public void inputParameters(ActionEvent event) {
        if (inputDataInterface == null) {
            inputDataInterface = new InputDataInterface(this);
        }
        inputDataInterface.inputParameters();
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
            xychart.setBounds(0.0, 1.25, 0.0, 200.0, 0.25, 50.0);
            xLowerBoundTextField.setText("0.0");
            xUpperBoundTextField.setText("1.25");
            yLowerBoundTextField.setText("0.0");
            yUpperBoundTextField.setText("200.0");
            xTickTextField.setText("0.25");
            yTickTextField.setText("50.0");
        } else if ((simControls instanceof CESTControls)) {
            xychart.setNames("CEST", "Offset (Hz)", "I(t)/I(0)", "20");
            xychart.setBounds(-6000, 6000, 0.0, 1.0, 2000.0, 0.25);
            xLowerBoundTextField.setText("-6000.0");
            xUpperBoundTextField.setText("6000.0");
            if (currentResProps != null) {
                double[] xVals = currentResProps.getExperimentData().stream().findFirst().get().getXVals();
                xychart.setBounds(Math.round(xVals[1]/1000)*1000, Math.round(xVals[xVals.length-1]/1000)*1000, 0.0, 1.0, 2000.0, 0.25);
                xLowerBoundTextField.setText(String.valueOf(Math.round(xVals[1]/1000)*1000));
                xUpperBoundTextField.setText(String.valueOf(Math.round(xVals[xVals.length-1]/1000)*1000));
            }
            yLowerBoundTextField.setText("0.0");
            yUpperBoundTextField.setText("1.0");
            xTickTextField.setText("2000.0");
            yTickTextField.setText("0.25");
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
        xychart.autoscaleBounds();

    }

    public void updateChartEquations(String equationName, double[] pars, double[] errs, double[] fields) {
        List<PlotEquation> equations = new ArrayList<>();
        for (int i = 0; i < fields.length; i++) {
            double[] extras = {fields[i] / fields[0]};
            //System.out.println("updateChartEquations got called with extras length = "+extras.length);
            PlotEquation plotEquation = new PlotEquation(equationName, pars, errs, extras);
            equations.add(plotEquation);
        }
        showEquations(equations);
    }

    public void showEquations(List<PlotEquation> equations) {
        xychart.setEquations(equations);
        xychart.layoutPlotChildren();
//        Optional<Double> rms = rms();

    }

    public XYBarChart getActiveChart() {
        return activeChart;
    }

    public XYBarChart addChart() {
        XYBarChart newChart = new XYBarChart();
        int nChildren = chartBox.getChildren().size();
        chartBox.getChildren().add(0, newChart);
        VBox.setVgrow(newChart, Priority.ALWAYS);
        activeChart = newChart;
        return newChart;

    }

    public void removeChart(Event e) {
        if ((activeChart != null) && (chartBox.getChildren().size() > 2)) {
            chartBox.getChildren().remove(activeChart);
            activeChart = (XYBarChart) chartBox.getChildren().get(0);
        }
    }

    public void updateTable(List<ResidueData> resDatas) {
        ObservableList<ResidueData.DataValue> data = FXCollections.observableArrayList();
        for (ResidueData resData : resDatas) {
            data.addAll(resData.getDataValues());
        }
        resInfoTable.itemsProperty().setValue(data);

        if (getFittingMode().equals("cpmg")) {
            TableColumn<ResidueData.DataValue, String> nameColumn = new TableColumn<>("Name");
            TableColumn<ResidueData.DataValue, String> resColumn = new TableColumn<>("Residue");
            TableColumn<ResidueData.DataValue, Double> xColumn = new TableColumn<>("Vcpmg");
            TableColumn<ResidueData.DataValue, Double> yColumn = new TableColumn<>("Reff");
            TableColumn<ResidueData.DataValue, String> errColumn = new TableColumn<>("Error");
            TableColumn<ResidueData.DataValue, String> peakColumn = new TableColumn<>("Peak");

            nameColumn.setCellValueFactory(new PropertyValueFactory<>("Name"));
            resColumn.setCellValueFactory(new PropertyValueFactory<>("Residue"));
            xColumn.setCellValueFactory(new PropertyValueFactory<>("X0"));
            yColumn.setCellValueFactory(new PropertyValueFactory<>("Y"));
            errColumn.setCellValueFactory(new PropertyValueFactory<>("Error"));
            peakColumn.setCellValueFactory(new PropertyValueFactory<>("Peak"));

            resInfoTable.getColumns().clear();
            resInfoTable.getColumns().addAll(nameColumn, resColumn, xColumn, yColumn, errColumn, peakColumn);
        } else if (getFittingMode().equals("cest")) {
            TableColumn<ResidueData.DataValue, String> nameColumn = new TableColumn<>("Name");
            TableColumn<ResidueData.DataValue, String> resColumn = new TableColumn<>("Residue");
            TableColumn<ResidueData.DataValue, Double> x0Column = new TableColumn<>("Offset");
            TableColumn<ResidueData.DataValue, Double> x1Column = new TableColumn<>("B1 Field");
            TableColumn<ResidueData.DataValue, Double> yColumn = new TableColumn<>("Intensity");
            TableColumn<ResidueData.DataValue, String> errColumn = new TableColumn<>("Error");

            nameColumn.setCellValueFactory(new PropertyValueFactory<>("Name"));
            resColumn.setCellValueFactory(new PropertyValueFactory<>("Residue"));
            x0Column.setCellValueFactory(new PropertyValueFactory<>("X0"));
            x1Column.setCellValueFactory(new PropertyValueFactory<>("X1"));
            yColumn.setCellValueFactory(new PropertyValueFactory<>("Y"));
            errColumn.setCellValueFactory(new PropertyValueFactory<>("Error"));

            resInfoTable.getColumns().clear();
            resInfoTable.getColumns().addAll(nameColumn, resColumn, x0Column, x1Column, yColumn, errColumn);
        }
    }

    public void updateTableWithPars(String mapName, String[] residues, String equationName, String state, List<int[]> allStates) {
        updateTableWithPars(mapName, residues, equationName, state, allStates, true);
    }

    public void updateTableWithPars(String mapName, String[] residues, String equationName, String state, List<int[]> allStates, boolean savePars) {
        List<ParValueInterface> allParValues = new ArrayList<>();
        if (residues != null) {
            for (String resNum : residues) {
                ResidueInfo resInfo = ChartUtil.getResInfo(mapName, String.valueOf(resNum));
                if (resInfo != null) {
                    currentResInfo = resInfo;
                    final String useEquationName;
                    if (equationName.equals("best")) {
                        useEquationName = resInfo.getBestEquationName();
                    } else {
                        useEquationName = equationName;
                    }
                    List<ParValueInterface> parValues = resInfo.getParValues(useEquationName, state);
                    if (resNum.equals(residues[0])) {
                        simControls.updateStates(allStates);
                        simControls.updateSliders(parValues, useEquationName);
                    }

                    allParValues.addAll(parValues);
                    CurveFit curveSet = currentResInfo.getCurveSet(useEquationName, state.replace("*", "0"));
//                    System.out.println("curv " + useEquationName + " " + state + " " + curveSet);
                    try {
                        String aic = String.format("%.2f", curveSet.getParMap().get("AIC"));
                        String rms = String.format("%.3f", curveSet.getParMap().get("RMS"));
                        aicLabel.setText(aic);
                        rmsLabel.setText(rms);
                    } catch (NullPointerException npEaic) {

                    }
                }
            }
            if (savePars) {
                currentMapName = mapName;
                currentResidues = new String[residues.length];
                System.arraycopy(residues, 0, currentResidues, 0, residues.length);
                currentState = state;
                currentEquationName = equationName;
                currentStates.clear();
                currentStates.addAll(allStates);
            }

            updateTableWithPars(allParValues);
        }
    }

    public void updateEquation(String mapName, String[] residues, String equationName) {
        String setEquation = "";
        for (String resNum : residues) {
            final String useEquationName;
            ResidueInfo resInfo = ChartUtil.getResInfo(mapName, resNum);
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
        equationChoice.setUserData("disabled");
        equationChoice.setValue(setEquation);
        equationChoice.setUserData(null);
    }

    public void updateTableWithPars(List<ParValueInterface> parValues) {
        DecimalFormat df = new DecimalFormat();
        ObservableList<ParValueInterface> data = FXCollections.observableArrayList();
        data.addAll(parValues);
        parameterTable.itemsProperty().setValue(data);

        TableColumn<ParValueInterface, String> residueColumn = new TableColumn<>("Residue");
        TableColumn<ParValueInterface, String> stateColumn = new TableColumn<>("State");
        TableColumn<ParValueInterface, String> nameColumn = new TableColumn<>("Name");
        TableColumn<ParValueInterface, String> valueColumn = new TableColumn<>("Value");
        TableColumn<ParValueInterface, String> errorColumn = new TableColumn<>("Error");

        residueColumn.setCellValueFactory(new PropertyValueFactory<>("Residue"));
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
        parameterTable.getColumns().addAll(residueColumn, stateColumn, nameColumn, valueColumn, errorColumn);
        ArrayList indices = new ArrayList<>();
        for (int i=0; i<parValues.size(); i++) {
            String state = stateColumn.getCellObservableValue(parValues.get(i)).getValue();
            String parName = nameColumn.getCellObservableValue(parValues.get(i)).getValue();
            if (state != null & parName != null) {
                if (!state.equals("0:0:0") & parName.equals("Kex") || !state.equals("0:0:0") & parName.equals("pA")) {
                    indices.add(i);
                }
            }
        }
        for (int i=0; i<indices.size(); i++) {
            parameterTable.getItems().remove((int) indices.get(i)-i);
        }
    }

    public void selectTableRow(String seriesName, int index) {
        parTabPane.getSelectionModel().select(1);
        List<ResidueData.DataValue> data = resInfoTable.getItems();
        int iRow = 0;
        for (ResidueData.DataValue dValue : data) {
            if ((dValue.getIndex() == index) && ((dValue.getName() + ":" + dValue.getResidue()).equals(seriesName))) {
                resInfoTable.getSelectionModel().clearAndSelect(iRow);
                resInfoTable.scrollTo(iRow);
                return;
            }
            iRow++;

        }

    }

    public void clearChart(Event e) {
        clearChart();
    }

    public void clearChart() {
        activeChart.getData().clear();
    }

    public void setYAxisType(String setName, String eqnName, String state, String typeName) {
        ObservableList<XYChart.Series<Double, Double>> data = ChartUtil.getParMapData(setName, eqnName, state, typeName);
        XYBarChart chart = activeChart;

        ((XYChart) chart).getData().addAll(data);
        if (chart.currentSeriesName.equals("")) {
            chart.currentSeriesName = data.get(0).getName();
        }
        NumberAxis xAxis = (NumberAxis) chart.getXAxis();
        NumberAxis yAxis = (NumberAxis) chart.getYAxis();
        xAxis.setAutoRanging(false);
        double xMin = Math.floor((ChartUtil.minRes - 2) / 5.0) * 5.0;
        double xMax = Math.ceil((ChartUtil.maxRes + 2) / 5.0) * 5.0;
        xAxis.setLowerBound(xMin);
        xAxis.setUpperBound(xMax);
        xAxis.setTickUnit(10);
        xAxis.setMinorTickCount(5);
        yAxis.setMinWidth(75.0);
        yAxis.setPrefWidth(75.0);

        if (typeName.equals("Kex")) {
            chart.getYAxis().setLabel("Kex");
            chart.setUnifyYAxes(true);
        } else {
            chart.getYAxis().setLabel(typeName);
            chart.setUnifyYAxes(false);
        }
        chart.setTitle(setName);
        if (chart.getData().size() > 1) {
            chart.setLegendVisible(true);
            chart.setLegendSide(Side.TOP);
        } else {
            chart.setLegendVisible(false);
        }
        currentResProps = ChartUtil.residueProperties.get(setName);
        chart.setResProps(currentResProps);
        if (currentResProps != null) {
            defaultField = currentResProps.getFields()[0];
        }
    }

    void makeAxisMenu() {
        axisMenu.getItems().clear();
        //Residue	 Peak	GrpSz	Group	Equation	   RMS	   AIC	Best	     R2	  R2.sd	    Rex	 Rex.sd	    Kex	 Kex.sd	     pA	  pA.sd	     dW	  dW.sd

        String[] parTypes = getParTypes();
        Map<String, ResidueProperties> residueProps = ChartUtil.residueProperties;
        // MenuItem clearItem = new MenuItem("Clear");
        // clearItem.setOnAction(e -> clearChart());
        // axisMenu.getItems().add(clearItem);
        for (String parType : parTypes) {
            Menu cascade = new Menu(parType);
            axisMenu.getItems().add(cascade);
            for (ResidueProperties residueProp : residueProps.values()) {
                String setName = residueProp.getName();
                Menu cascade2 = new Menu(setName);
                cascade.getItems().add(cascade2);
                ArrayList<String> equationNames = new ArrayList<>();
                equationNames.add("best");
                equationNames.addAll(residueProp.getEquationNames());
                List<String> stateStrings = residueProp.getStateStrings();

                for (String equationName : equationNames) {
                    if ((stateStrings.size() < 2) || parType.equals("RMS") || parType.equals("AIC") || parType.equals("Equation")) {
                        MenuItem cmItem1 = new MenuItem(equationName);
                        cmItem1.setOnAction(e -> setYAxisType(setName, equationName, "0:0:0", parType));
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
                                cmItem1.setOnAction(e -> setYAxisType(setName, equationName, state, parType));
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
        if (currentResProps == null) {
            guessSimData();
        } else {
            try {
                EquationFitter equationFitter = getFitter();
                String[] resNums = {String.valueOf(currentResInfo.getResNum())};
                equationFitter.setData(currentResProps, resNums);
                String equationName = simControls.getEquation();
//        System.out.println("guesses eqnFitter = " + equationFitter);
//        System.out.println("guesses resNums = " + resNums);
//        System.out.println("guesses eqnName = " + equationName);
                List<ParValueInterface> guesses = equationFitter.guessPars(equationName, absValueModeCheckBox.isSelected());
                simControls.updateSliders(guesses, equationName);
                simControls.simSliderAction("");
            } catch (NullPointerException npE1) {
                Alert alert = new Alert(Alert.AlertType.ERROR);
                alert.setContentText("Error: Residue must be selected");
                alert.showAndWait();
                return;
            }
        }
    }

    public Optional<Double> rms() {
        Optional<Double> rms;
        try {
            EquationFitter equationFitter = getFitter();
            if (currentResInfo != null) {
                String[] resNums = {String.valueOf(currentResInfo.getResNum())};
                equationFitter.setData(currentResProps, resNums);
                String equationName = simControls.getEquation();
                equationFitter.setupFit(equationName, absValueModeCheckBox.isSelected());
                int[][] map = equationFitter.getFitModel().getMap();
                double[] sliderGuesses = simControls.sliderGuess(equationName, map);
                rms = Optional.of(equationFitter.rms(sliderGuesses));
            } else {
                rms = Optional.empty();
            }
        } catch (NullPointerException npE2) {
            rms = Optional.empty();
        }
        return rms;
    }

    @FXML
    public void fitEquation(ActionEvent event) {
//        EquationFitter equationFitter = new CPMGFit();
        fitResult = null;
        try {
            EquationFitter equationFitter = getFitter();
            if (currentResProps == null) {
                fitSimData();
            } else {
                String[] resNums = {String.valueOf(currentResInfo.getResNum())};
                equationFitter.setData(currentResProps, resNums);
                String equationName = simControls.getEquation();
                equationFitter.setupFit(equationName, absValueModeCheckBox.isSelected());
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
                fitResult = equationFitter.doFit(equationName, absValueModeCheckBox.isSelected(), nonParBootStrapCheckBox.isSelected(), sliderGuesses);
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

    public void updateAfterFit(CPMGFitResult fitResult) {
        List<PlotEquation> equations = new ArrayList<>();
        int nCurves = fitResult.getNCurves();
        for (int iCurve = 0; iCurve < nCurves; iCurve++) {
            CurveFit curveFit = fitResult.getCurveFit(iCurve);
            List<ParValueInterface> parValues = curveFit.getParValues();
            equationChoice.getSelectionModel().select(fitResult.getEquationName());
            String aic = String.format("%.2f", fitResult.getAicc());
            String rms = String.format("%.3f", fitResult.getRms());
            aicLabel.setText(aic);
            rmsLabel.setText(rms);
            updateTableWithPars(parValues);
            simControls.updateSliders(parValues, fitResult.getEquationName());
            String equationName = fitResult.getEquationName(); //equationSelector.getValue();
            //System.out.println("Fit button residueProperties = " + residueProperties);
            //System.out.println("Fit button expData = " + residueProps.getExperimentData("cest"));
            Optional<ExperimentData> optionalData = Optional.empty();
            if (currentResProps != null) {
                optionalData = currentResProps.getExperimentData().stream().findFirst();
            }

            if (optionalData.isPresent() && optionalData.get().getExtras().size() > 0) {
                ExperimentData expData = optionalData.get();
                double[] pars = curveFit.getEquation().getPars(); //pars = getPars(equationName);
                double[] errs = curveFit.getEquation().getErrs(); //double[] errs = new double[pars.length];
                double[] extras = new double[3];
                for (int j = 0; j < expData.getExtras().size() / 2; j++) {
                    extras[0] = expData.getField();
                    extras[1] = expData.getExtras().get(2 * j);
                    extras[2] = expData.getExtras().get(2 * j + 1);
//                    System.out.println("Fit button expData extras size = " + expData.getExtras().size() + " extra[1] = " + extras[1]);
                    PlotEquation plotEquation = new PlotEquation(equationName, pars, errs, extras);
                    //equationCopy.setExtra(extras);

                    equations.add(plotEquation);
                }
            } else {
                double[] pars = curveFit.getEquation().getPars(); //pars = getPars(equationName);
                double[] errs = curveFit.getEquation().getErrs(); //double[] errs = new double[pars.length];
                double[] extras = new double[3];
                double[] simExtras = simControls.getExtras();
                extras[0] = CPMGFit.REF_FIELD;

                for (int i = 0; i < simExtras.length; i++) {
                    extras[i + 1] = simExtras[i];
                }
                PlotEquation plotEquation = new PlotEquation(equationName, pars, errs, extras);

                //equationCopy.setExtra(extras);
                //System.out.println("Fit button expData extras size = " + expData.getExtras().size() + " extra[0] = " + extras[0]);
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
        residueFitter.fitResidues(currentResProps);
//        }
    }

    @FXML
    public void fitGroupResidues(ActionEvent event) {
//        if (getFittingMode().equals("cest")) {
//            ChooseCESTFitEquations.allRes = false;
//            ChooseCESTFitEquations.create();
//        } else {
        fitResult = null;
        List<List<String>> allResidues = new ArrayList<>();
        List<String> groupResidues = new ArrayList<>();
        XYBarChart chart = getActiveChart();
        groupResidues.addAll(chart.selectedResidues);
        if (!groupResidues.isEmpty()) {
            allResidues.add(groupResidues);
            currentResProps.setAbsValueMode(absValueModeCheckBox.isSelected());
            if (nonParBootStrapCheckBox.isSelected()) {
                currentResProps.setBootStrapMode("nonparametric");
            } else {
                currentResProps.setBootStrapMode("parametric");
            }
            residueFitter.fitResidues(currentResProps, allResidues);
        }
//        }
    }

    public void refreshFit() {
        XYBarChart chart = getActiveChart();
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
        File file = fileChooser.showSaveDialog(MainApp.primaryStage);
        if (file != null) {
            DataIO.saveParametersToFile(file.getAbsolutePath(), currentResProps);
        }
    }

    public Double updateFitProgress(Double f) {
        XYBarChart chart = getActiveChart();
        String seriesName = chart.currentSeriesName;
        String[] sParts;
        if (seriesName.length() == 0) {
            sParts = new String[4];
            sParts[0] = currentResProps.getName();
            sParts[1] = "best";
            sParts[2] = "0:0:0";
            sParts[3] = "RMS";
        } else {
            sParts = seriesName.split("\\|");
        }
        if (Platform.isFxApplicationThread()) {
            clearChart();
            statusBar.setProgress(f);
            setYAxisType(sParts[0], sParts[1], sParts[2], sParts[3]);
            currentResProps = ChartUtil.residueProperties.get(currentResProps.getName());

        } else {
            Platform.runLater(() -> {
                clearChart();
                setYAxisType(sParts[0], sParts[1], sParts[2], sParts[3]);
                statusBar.setProgress(f);
                currentResProps = ChartUtil.residueProperties.get(currentResProps.getName());

            });
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
            snapit(xychart, file);
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
                    if (node instanceof XYBarChart) {
                        XYBarChart barChart = (XYBarChart) node;
                        String barChartName = barChart.toString();
                        HashMap<String, Object> chartData = barChart.getChartData();
                        barChartData.add(chartData);
                    }
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

    public CPMGFitResult getFitResult() {
        return fitResult;
    }

    public ResidueInfo getResidueInfo() {
        return currentResInfo;
    }

    public String getEquationName() {
        return currentEquationName;
    }

    public String getFittingMode() {
        String fitMode = "cpmg";
        if (currentResProps == null) {
            String simMode = getSimMode();
            if (simMode != null) {
                fitMode = simMode;
            }
        } else {
            fitMode = currentResProps.getExpMode();
        }
        return fitMode.toLowerCase();
    }

    public EquationFitter getFitter() {
        if (getFittingMode().equals("exp")) {
            return new ExpFit();
        } else if (getFittingMode().equals("cpmg")) {
            return new CPMGFit();
        } else if (getFittingMode().equals("cest")) {
            return new CESTFit();
        }
        return null;
    }

    public String[] getParNames(String equationName) {
        EquationType type = ResidueFitter.getEquationType(equationName);
        return type.getParNames();
    }

    public String[] getParTypes() {
        String[] cpmgTypes = {"R2", "Rex", "Kex", "pA", "dW", "RMS", "AIC", "Equation"};
        String[] expTypes = {"A", "R", "C", "RMS", "AIC", "Equation"};
        String[] cestTypes = {"kex", "pb", "deltaA0", "deltaB0", "R1A", "R1B", "R2A", "R2B", "RMS", "AIC", "Equation"};
        String[] nullTypes = {"RMS", "AIC", "Equation"};
        if (getFittingMode().equals("exp")) {
            return expTypes;
        } else if (getFittingMode().equals("cpmg")) {
            return cpmgTypes;
        } else if (getFittingMode().equals("cest")) {
            return cestTypes;
        }
        return nullTypes;
    }

    @FXML
    void setBestEquation(ActionEvent e) {
        for (String resNum : currentResidues) {
            ResidueInfo resInfo = ChartUtil.getResInfo(currentMapName, String.valueOf(resNum));
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
            if (!currentStates.isEmpty() && equationName != null) {
                // copy it so it doesn't get cleared by clear call in updateTableWithPars
                List<int[]> useStates = new ArrayList<>(currentStates);
                updateTableWithPars(currentMapName, currentResidues, equationName, currentState, useStates, false);
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
        currentResidues = null;
        currentResInfo = null;
        currentResProps = null;
        fitResult = null;
        xychart.clear();
        chartBox.getChildren().remove(activeChart);
        activeChart = xyBarChart0;
        getSimMode();
        setSimControls();
        updateXYChartLabels();
        simControls.simSliderAction("");
    }

    void showInfo(String equationName) {
        String mapName = currentMapName;
        String state = currentState;
        showInfo(currentResProps, equationName, mapName, state, currentResidues, xychart);
    }

    void showInfo(ResidueProperties resProps, String equationName, String mapName, String state, String[] residues, PlotData plotData) {
        ArrayList<PlotEquation> equations = new ArrayList<>();
        ObservableList<XYChart.Series<Double, Double>> allData = FXCollections.observableArrayList();
        List<ResidueData> resDatas = new ArrayList<>();
        List<int[]> allStates = new ArrayList<>();
        if (resProps != null) {
            for (ExperimentData expData : resProps.getExperimentData()) {
                if (residues == null) {
                    break;
                }
                for (String residue : residues) {
//                    System.out.println("get resd " + residue + " " + expData.getResidueData(residue));

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
        }
        plotData.setData(allData);
        setBounds();
        plotData.setEquations(equations);
        plotData.layoutPlotChildren();
        updateTable(resDatas);
        if (residues != null) {
            updateTableWithPars(mapName, residues, equationName, state, allStates);
            updateEquation(mapName, residues, equationName);
        }

    }

    public void setSimData(EquationFitter equationFitter) {
        ArrayList<Double> xValues = getSimXData();
        ArrayList<Double> yValues = getSimYData();
        double[] extras = simControls.getExtras();
        ArrayList[] allXValues = new ArrayList[extras.length + 1];
        allXValues[0] = xValues;
        ArrayList<Double> errValues = new ArrayList<>();
        for (int i = 0; i < yValues.size(); i++) {
            errValues.add(yValues.get(i) * 0.05);
        }
        for (int j = 0; j < extras.length; j++) {
            ArrayList<Double> xValuesEx = new ArrayList<>();
            for (int i = 0; i < yValues.size(); i++) {
                xValuesEx.add(extras[j]);
            }
            allXValues[1 + j] = xValuesEx;
        }
        equationFitter.setData(allXValues, yValues, errValues);
    }

    public void guessSimData() {
        EquationFitter equationFitter = getFitter();
        setSimData(equationFitter);
        String equationName = simControls.getEquation();
        List<ParValueInterface> guesses = equationFitter.guessPars(equationName, absValueModeCheckBox.isSelected());
        simControls.updateSliders(guesses, equationName);
        simControls.simSliderAction("");
    }

    public void fitSimData() {
        EquationFitter equationFitter = getFitter();
        setSimData(equationFitter);
        String equationName = simControls.getEquation();
        equationFitter.setupFit(equationName, absValueModeCheckBox.isSelected());
        int[][] map = equationFitter.getFitModel().getMap();
        double[] sliderGuesses = null;
        if (sliderGuessCheckBox.isSelected()) {
            sliderGuesses = simControls.sliderGuess(equationName, map);
        }
        fitResult = equationFitter.doFit(equationName, absValueModeCheckBox.isSelected(), nonParBootStrapCheckBox.isSelected(), sliderGuesses);
        updateAfterFit(fitResult);
    }

    ArrayList<Double> getSimXData() {
        ObservableList<XYChart.Series<Double, Double>> allData = xychart.getData();
        ArrayList<Double> xValues = new ArrayList<>();
        for (XYChart.Series<Double, Double> series : allData) {
            for (XYChart.Data<Double, Double> dataPoint : series.getData()) {
                double x = dataPoint.getXValue();
                xValues.add(x);
            }
        }
        return xValues;
    }

    ArrayList<Double> getSimYData() {
        ObservableList<XYChart.Series<Double, Double>> allData = xychart.getData();
        ArrayList<Double> yValues = new ArrayList<>();
        for (XYChart.Series<Double, Double> series : allData) {
            for (XYChart.Data<Double, Double> dataPoint : series.getData()) {
                double y = dataPoint.getYValue();
                yValues.add(y);
            }
        }
        return yValues;
    }

    @FXML
    void showSimData(ActionEvent e) {
        ObservableList<XYChart.Series<Double, Double>> allData = FXCollections.observableArrayList();
        String equationName = simControls.getEquation();
        EquationType eType = ResidueFitter.getEquationType(equationName);
        int[][] map = eType.makeMap(1);
        EquationFitter equationFitter = getFitter();
        equationFitter.getFitModel().setMap(map);
        double[] sliderGuesses = simControls.sliderGuess(equationName, map);
        double[] xBounds = xychart.getXBounds();
        double[] yBounds = xychart.getYBounds();
        double sdev = Math.abs(yBounds[1] - yBounds[0]) * 0.02;
        double[] xValues = equationFitter.getSimX();
        double fieldRef = 1.0;
        int iLine = 0;

        List<XYChart.Series<Double, Double>> data = new ArrayList<>();
        XYChart.Series<Double, Double> series = new XYChart.Series<>();
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
                double yValue = eqn.calculate(sliderGuesses, ax, eqn.getExtra(0));
                yValue += sdev * rand.nextGaussian();
                XYChart.Data dataPoint = new XYChart.Data(xValue, yValue);
                dataPoint.setExtraValue(new Double(sdev));
                series.getData().add(dataPoint);
            }
        }
        allData.addAll(data);
        xychart.setData(allData);
    }

    @FXML
    public void loadSimData(ActionEvent e) {
        ObservableList<XYChart.Series<Double, Double>> allData = FXCollections.observableArrayList();
        FileChooser fileChooser = new FileChooser();
        File file = fileChooser.showOpenDialog(null);
        if (file != null) {
            try {
                List<Double>[] dataValues = DataIO.loadSimData(file);
                List<XYChart.Series<Double, Double>> data = new ArrayList<>();
                XYChart.Series<Double, Double> series = new XYChart.Series<>();
                series.setName("sim" + ":" + "0");
                data.add(series);
                for (int i = 0; i < dataValues[0].size(); i++) {
                    XYChart.Data dataPoint = new XYChart.Data(dataValues[0].get(i), dataValues[1].get(i));
                    if (!dataValues[2].isEmpty()) {
                        dataPoint.setExtraValue(new Double(dataValues[2].get(i)));
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
        Platform.exit();
        System.exit(0);
    }

}

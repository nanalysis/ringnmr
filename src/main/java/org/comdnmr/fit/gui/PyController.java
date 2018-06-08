package org.comdnmr.fit.gui;

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
import java.io.PrintStream;
import java.util.HashMap;
import javafx.scene.control.Alert;
import javafx.scene.control.TextArea;
import javafx.scene.control.SplitPane;
import static org.comdnmr.fit.gui.ChartUtil.residueProperties;
import javafx.scene.Scene;
import javafx.scene.chart.ScatterChart;
import javafx.scene.layout.HBox;
import javafx.scene.control.Label;
import javafx.scene.control.TextField;
import javafx.scene.control.Button;

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
    TextArea textArea;

    @FXML
    ChoiceBox<String> simChoice;

    @FXML
    CheckBox sliderGuessCheckBox;

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

    VBox MCdisplay = new VBox();
    Scene stageScene = new Scene(MCdisplay, 500, 500);
    ChoiceBox<String> MCxArrayChoice = new ChoiceBox<>();
    ChoiceBox<String> MCyArrayChoice = new ChoiceBox<>();
    XYChart activeMCchart;
    XYChart.Series MCseries = new XYChart.Series();

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
            alert.setContentText("Error: Residue must be selected");
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

        PrintStream printStream = new PrintStream(new ConsoleRedirect(textArea));
        // keeps reference of standard output stream
        //standardOut = System.out;
        //standardErr = System.err;

        // re-assigns standard output stream and error output stream
        System.setOut(printStream);
        System.setErr(printStream);

        simChoice.getItems().add("Simulate CPMG");
        simChoice.getItems().add("Simulate EXP");
        simChoice.getItems().add("Simulate CEST");
        simChoice.setValue("Simulate CPMG");
        simChoice.valueProperty().addListener(s -> {
            simAction();
        });

        splitPane.setDividerPositions(0.4, 0.7);

        setBoundsButton.setOnAction(this::setBounds);
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

    public void clearConsole() {
        textArea.setText("");
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
            xychart.setBounds(-6000.0, 6000.0, 0.0, 1.0, 2000.0, 0.25);
            xLowerBoundTextField.setText("-6000.0");
            xUpperBoundTextField.setText("6000.0");
            yLowerBoundTextField.setText("0.0");
            yUpperBoundTextField.setText("1.0");
            xTickTextField.setText("2000.0");
            yTickTextField.setText("0.25");
        }
    }

    public void setBounds(ActionEvent event) {
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
        ObservableList<ParValueInterface> data = FXCollections.observableArrayList();
        data.addAll(parValues);
        parameterTable.itemsProperty().setValue(data);

        TableColumn<ParValueInterface, String> residueColumn = new TableColumn<>("Residue");
        TableColumn<ParValueInterface, String> stateColumn = new TableColumn<>("State");
        TableColumn<ParValueInterface, String> nameColumn = new TableColumn<>("Name");
        TableColumn<ParValueInterface, Double> valueColumn = new TableColumn<>("Value");
        TableColumn<ParValueInterface, Double> errorColumn = new TableColumn<>("Error");

        residueColumn.setCellValueFactory(new PropertyValueFactory<>("Residue"));
        stateColumn.setCellValueFactory(new PropertyValueFactory<>("State"));
        nameColumn.setCellValueFactory(new PropertyValueFactory<>("Name"));
        valueColumn.setCellValueFactory(new PropertyValueFactory<>("Value"));
        errorColumn.setCellValueFactory(new PropertyValueFactory<>("Error"));

        parameterTable.getColumns().clear();
        parameterTable.getColumns().addAll(residueColumn, stateColumn, nameColumn, valueColumn, errorColumn);
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
        try {
            EquationFitter equationFitter = getFitter();
            String[] resNums = {String.valueOf(currentResInfo.getResNum())};
            equationFitter.setData(currentResProps, resNums);
            String equationName = simControls.getEquation();
//        System.out.println("guesses eqnFitter = " + equationFitter);
//        System.out.println("guesses resNums = " + resNums);
//        System.out.println("guesses eqnName = " + equationName);
            List<ParValueInterface> guesses = equationFitter.setupFit(equationName, absValueModeCheckBox.isSelected());
            simControls.updateSliders(guesses, equationName);
            simControls.simSliderAction("");
        } catch (NullPointerException npE1) {
            Alert alert = new Alert(Alert.AlertType.ERROR);
            alert.setContentText("Error: Residue must be selected");
            alert.showAndWait();
            return;
        }
    }

    @FXML
    public void fitEquation(ActionEvent event) {
//        EquationFitter equationFitter = new CPMGFit();
        fitResult = null;
        try {
            EquationFitter equationFitter = getFitter();
            String[] resNums = {String.valueOf(currentResInfo.getResNum())};
            equationFitter.setData(currentResProps, resNums);
            String equationName = simControls.getEquation();
//            equationFitter.getFitModel().setMap(stateCount, states);
//            int[][] map = equationFitter.getFitModel().getMap();
//            System.out.println("map = " + map);
//            System.out.println("getFitModel = " + equationFitter.getFitModel());
//        System.out.println("fitEqn eqnFitter = " + equationFitter);
//        System.out.println("fitEqn resNums = " + resNums);
//        System.out.println("fitEqn eqnName = " + equationName);
            double[] sliderGuesses = null; //simControls.sliderGuess(equationName, map);
            fitResult = equationFitter.doFit(equationName, absValueModeCheckBox.isSelected(), nonParBootStrapCheckBox.isSelected(), sliderGuesses);
            updateAfterFit(fitResult);
        } catch (NullPointerException npE2) {
            Alert alert = new Alert(Alert.AlertType.ERROR);
            alert.setContentText("Error: Residue must be selected");
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
            updateTableWithPars(parValues);
            simControls.updateSliders(parValues, fitResult.getEquationName());
            String equationName = fitResult.getEquationName(); //equationSelector.getValue();
            //System.out.println("Fit button residueProperties = " + residueProperties);
            ResidueProperties residueProps = residueProperties.get("cest"); // fixme
            //System.out.println("Fit button expData = " + residueProps.getExperimentData("cest"));
            ExperimentData expData = null;
            if (residueProps != null) {
                expData = residueProps.getExperimentData("cest"); // fixme
            }

            if (expData != null && expData.getExtras().size() > 0) {
                double[] pars = curveFit.getEquation().getPars(); //pars = getPars(equationName);
                double[] errs = curveFit.getEquation().getErrs(); //double[] errs = new double[pars.length];
                double[] extras = new double[2];
                for (int j = 0; j < expData.getExtras().size(); j++) {
                    extras[0] = 1.0;
                    extras[1] = expData.getExtras().get(j) * 2 * Math.PI;
                    //System.out.println("Fit button expData extras size = " + expData.getExtras().size() + " extra[1] = " + extras[1]);
                    PlotEquation plotEquation = new PlotEquation(equationName, pars, errs, extras);
                    //equationCopy.setExtra(extras);

                    equations.add(plotEquation);
                }
            } else {
                double[] pars = curveFit.getEquation().getPars(); //pars = getPars(equationName);
                double[] errs = curveFit.getEquation().getErrs(); //double[] errs = new double[pars.length];
                double[] extras = new double[1];
                extras[0] = 1.0;
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
        fitResult = null;
        residueFitter.fitResidues(currentResProps);
    }

    @FXML
    public void fitGroupResidues(ActionEvent event) {
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

    public void showMCplot(ActionEvent event) {
        //Create new Stage for popup window  
        Stage MCStage = new Stage();
        MCStage.setTitle("Monte Carlo Results: " + simControls.getEquation());
        Label xlabel = new Label("  X Array:  ");
        Label ylabel = new Label("  Y Array:  ");

        //Populate ChoiceBoxes with fitting variable names
        MCxArrayChoice.getItems().clear();
        MCyArrayChoice.getItems().clear();
        MCxArrayChoice.getItems().addAll(simControls.getParNames());
        MCyArrayChoice.getItems().addAll(simControls.getParNames());
        MCxArrayChoice.setValue(simControls.getParNames().get(1));
        MCyArrayChoice.setValue(simControls.getParNames().get(0));

        try {
            MCxArrayChoice.valueProperty().addListener(x -> {
                updateMCplot();
            });
            MCyArrayChoice.valueProperty().addListener(y -> {
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
        hBox.getChildren().addAll(xlabel, MCxArrayChoice, ylabel, MCyArrayChoice);

        //Create the Scatter chart  
        NumberAxis MCxAxis = new NumberAxis();
        NumberAxis MCyAxis = new NumberAxis();
        ScatterChart<Number, Number> MCchart = new ScatterChart(MCxAxis, MCyAxis);
        MCxAxis.setAutoRanging(true);
        MCxAxis.setForceZeroInRange(false);
        MCyAxis.setAutoRanging(true);
        MCyAxis.setForceZeroInRange(false);

        MCxAxis.setLabel(MCxArrayChoice.getValue());
        MCyAxis.setLabel(MCyArrayChoice.getValue());

        //Prepare XYChart.Series objects by setting data 
        MCseries.getData().clear();
        int xInd = MCxArrayChoice.getSelectionModel().getSelectedIndex();
        int yInd = MCyArrayChoice.getSelectionModel().getSelectedIndex();
        if (fitResult != null && fitResult.getSimPars() != null) {
            for (int i = 0; i < fitResult.getSimPars()[xInd].length; i++) {
                MCseries.getData().add(new XYChart.Data(fitResult.getSimPars()[xInd][i], fitResult.getSimPars()[yInd][i]));
            }

            //Setting the data to scatter chart   
            MCchart.getData().clear();
            MCchart.getData().addAll(MCseries);

            activeMCchart = MCchart;

            MCdisplay.getChildren().clear();
            MCdisplay.getChildren().add(hBox);
            MCdisplay.getChildren().add(activeMCchart);

            MCStage.setScene(stageScene);
            MCStage.show();
        } else if (residueFitter.getFitResult() != null && residueFitter.getFitResult().getSimPars() != null) {
            for (int i = 0; i < residueFitter.getFitResult().getSimPars()[xInd].length; i++) {
                MCseries.getData().add(new XYChart.Data(residueFitter.getFitResult().getSimPars()[xInd][i], residueFitter.getFitResult().getSimPars()[yInd][i]));
            }

            //Setting the data to scatter chart   
            MCchart.getData().clear();
            MCchart.getData().addAll(MCseries);

            activeMCchart = MCchart;

            MCdisplay.getChildren().clear();
            MCdisplay.getChildren().add(hBox);
            MCdisplay.getChildren().add(activeMCchart);

            MCStage.setScene(stageScene);
            MCStage.show();
        } else {
            Alert alert = new Alert(Alert.AlertType.ERROR);
            alert.setContentText("Error: Fit must first be performed.");
            alert.showAndWait();
            return;
        }

    }

    public void updateMCplot() {
        try {
            //Create the Scatter chart 
            MCdisplay.getChildren().remove(activeMCchart);
            NumberAxis MCxAxis = new NumberAxis();
            NumberAxis MCyAxis = new NumberAxis();
            ScatterChart<Number, Number> MCchart = new ScatterChart(MCxAxis, MCyAxis);
            MCxAxis.setAutoRanging(true);
            MCxAxis.setForceZeroInRange(false);
            MCyAxis.setAutoRanging(true);
            MCyAxis.setForceZeroInRange(false);

            MCxAxis.setLabel(MCxArrayChoice.getSelectionModel().getSelectedItem());
            MCyAxis.setLabel(MCyArrayChoice.getSelectionModel().getSelectedItem());

            MCseries.getData().clear();
            int xInd = MCxArrayChoice.getSelectionModel().getSelectedIndex();
            int yInd = MCyArrayChoice.getSelectionModel().getSelectedIndex();
            if (fitResult != null && fitResult.getSimPars() != null) {
                double[][] simPars = fitResult.getSimPars();
                for (int i = 0; i < simPars[xInd].length; i++) {
                    MCseries.getData().add(new XYChart.Data(simPars[xInd][i], simPars[yInd][i]));
                }

                MCchart.getData().clear();
                MCchart.getData().addAll(MCseries);

                activeMCchart = MCchart;
                MCdisplay.getChildren().add(activeMCchart);
            } else if (residueFitter.getFitResult() != null && residueFitter.getFitResult().getSimPars() != null) {
                double[][] simPars = residueFitter.getFitResult().getSimPars();
                for (int i = 0; i < residueFitter.getFitResult().getSimPars()[xInd].length; i++) {
                    MCseries.getData().add(new XYChart.Data(simPars[xInd][i], simPars[yInd][i]));
                }

                MCchart.getData().clear();
                MCchart.getData().addAll(MCseries);

                activeMCchart = MCchart;
                MCdisplay.getChildren().add(activeMCchart);
            } else {
                Alert alert = new Alert(Alert.AlertType.ERROR);
                alert.setContentText("Error: Fit must first be performed.");
                alert.showAndWait();
                return;
            }
        } catch (ArrayIndexOutOfBoundsException aioobE) {
        }

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

    public void exportExecutable(String fileSuggestion, String exportType) throws ScriptException {
        FileChooser chooser = new FileChooser();
        chooser.setInitialFileName(fileSuggestion);
        File file = chooser.showSaveDialog(MainApp.primaryStage);
        if (file != null) {
            String filePath = file.getAbsolutePath();

            Map<String, Object>[] plottedData = xychart.getPlottedData();
            Map<String, Object> graphData = xychart.getGraphData();

            graphData.put("file", filePath);
            graphData.put("exportType", exportType);
            ArrayList<Object> barChartData;
            if (exportType != "grace") {
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
        exportExecutable("ASCII.agr", "grace");
    }

    public void savePythonFile() throws IOException, ScriptException {
        exportExecutable("graph.py", "python");
    }

    public void saveRFile() throws IOException, ScriptException {
        exportExecutable("graph.r", "r");
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

    public String getFittingMode() {
        if (currentResProps == null) {
            return "cpmg";
        } else {
            return currentResProps.getExpMode();
        }
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
            if (!currentStates.isEmpty()) {
                // copy it so it doesn't get cleared by clear call in updateTableWithPars
                List<int[]> useStates = new ArrayList<>(currentStates);
                updateTableWithPars(currentMapName, currentResidues, equationName, currentState, useStates, false);
                showInfo(equationName);
            }
        }
    }

    public String getSimMode() {
        String simSelected = simChoice.getSelectionModel().getSelectedItem();
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
                    System.out.println("get resd " + residue + " " + expData.getResidueData(residue));

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
        updateXYChartLabels();
        plotData.setEquations(equations);
        plotData.layoutPlotChildren();
        updateTable(resDatas);
        if (residues != null) {
            updateTableWithPars(mapName, residues, equationName, state, allStates);
            updateEquation(mapName, residues, equationName);
        }

    }

}

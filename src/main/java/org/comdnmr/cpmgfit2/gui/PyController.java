package org.comdnmr.cpmgfit2.gui;

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
import org.comdnmr.cpmgfit2.calc.CPMGFit;
import org.comdnmr.cpmgfit2.calc.CPMGFitResult;
import org.comdnmr.cpmgfit2.calc.CurveFit;
import org.comdnmr.cpmgfit2.calc.DataIO;
import org.comdnmr.cpmgfit2.calc.EquationFitter;
import org.comdnmr.cpmgfit2.calc.EquationType;
import org.comdnmr.cpmgfit2.calc.ExpFit;
import org.comdnmr.cpmgfit2.calc.ExperimentData;
import org.comdnmr.cpmgfit2.calc.ResidueData;
import org.controlsfx.control.PropertySheet;
import org.comdnmr.cpmgfit2.calc.ParValueInterface;
import org.comdnmr.cpmgfit2.calc.PlotEquation;
import org.comdnmr.cpmgfit2.calc.ProcessingStatus;
import org.comdnmr.cpmgfit2.calc.ResidueFitter;
import org.comdnmr.cpmgfit2.calc.ResidueInfo;
import org.comdnmr.cpmgfit2.calc.ResidueProperties;
import org.controlsfx.control.StatusBar;
import org.comdnmr.cpmgfit2.calc.CESTFit;
import java.io.PrintStream;
import javafx.scene.control.Alert;
import javafx.scene.control.TextArea;
import org.comdnmr.cpmgfit2.calc.ParValue;
import static org.comdnmr.cpmgfit2.gui.ChartUtil.residueProperties;

public class PyController implements Initializable {

    public static PyController mainController;
    XYBarChart activeChart;
    @FXML
    XYBarChart xyBarChart0;

    @FXML
    PlotData cpmgchart;

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

    EquationControls simControls;

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
        } else if (getFittingMode().equals("cest")) {
            simControls = new CESTControls();
            equationChoice.getItems().clear();
            equationChoice.getItems().add("+");
            equationChoice.getItems().addAll(CESTFit.getEquationNames());
            equationChoice.setValue(CESTFit.getEquationNames().get(0));
        } else {
            simControls = new ExpControls();
            equationChoice.getItems().clear();
            equationChoice.getItems().add("+");
            equationChoice.getItems().addAll(ExpFit.getEquationNames());
            equationChoice.setValue(ExpFit.getEquationNames().get(0));
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

    public void loadParameterFile(Event e) {
        FileChooser fileChooser = new FileChooser();
        fileChooser.getExtensionFilters().addAll(new FileChooser.ExtensionFilter("Yaml File", "*.yaml", "*.yml"));
        Stage stage = MainApp.primaryStage;
        File file = fileChooser.showOpenDialog(stage);
        if (file != null) {
            if (activeChart != null) {
                clearChart();
                currentResidues = null;
            }
            ChartUtil.loadParameters(file.toString());
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
        cpmgchart.setEquations(equations);
        cpmgchart.layoutPlotChildren();

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
        try {
            EquationFitter equationFitter = getFitter();
            String[] resNums = {String.valueOf(currentResInfo.getResNum())};
            equationFitter.setData(currentResProps, resNums);
            String equationName = simControls.getEquation();
//        System.out.println("fitEqn eqnFitter = " + equationFitter);
//        System.out.println("fitEqn resNums = " + resNums);
//        System.out.println("fitEqn eqnName = " + equationName);
            CPMGFitResult fitResult = equationFitter.doFit(equationName, absValueModeCheckBox.isSelected(), nonParBootStrapCheckBox.isSelected());
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
        residueFitter.fitResidues(currentResProps);
    }

    @FXML
    public void fitGroupResidues(ActionEvent event) {
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
            snapit(cpmgchart, file);
        }
    }

    public void saveGraceFile() throws IOException, ScriptException {
        FileChooser chooser = new FileChooser();
        chooser.setInitialFileName("ASCII.agr");
        File file = chooser.showSaveDialog(MainApp.primaryStage);
        if (file != null) {
            String filePath = file.getAbsolutePath();
            MainApp.engine.put("file", filePath);

            MainApp.engine.eval("writer = Writer()");
            int numPlots = cpmgchart.getNumPlots();
            int numFits = cpmgchart.getNumEquations();

            if (numPlots == numFits) {
                for (int i = 0; i < numPlots; i++) {
                    int color = i + 1;
                    MainApp.engine.put("color", color);

                    ArrayList<String> rawData = cpmgchart.returnLine(i);
                    MainApp.engine.put("rawData", rawData);
                    MainApp.engine.eval("writer.addRawData(rawData,color)");

                    ArrayList<String> fittedData = cpmgchart.returnEquation(i);
                    MainApp.engine.put("fittedData", fittedData);
                    MainApp.engine.eval("writer.addFittedData(fittedData,color)");
                }
            }
            MainApp.engine.eval("writer.write();");
            MainApp.engine.eval("writer.writeOut(file)");

            // TODO Figure out this part           
        }
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

    void showInfo(String equationName) {
        String mapName = currentMapName;
        String state = currentState;
        showInfo(currentResProps, equationName, mapName, state, currentResidues, cpmgchart);
    }

    void showInfo(ResidueProperties resProps, String equationName, String mapName, String state, String[] residues, PlotData plotData) {
        ArrayList<PlotEquation> equations = new ArrayList<>();
        ObservableList<XYChart.Series<Double, Double>> allData = FXCollections.observableArrayList();
        List<ResidueData> resDatas = new ArrayList<>();
        List<int[]> allStates = new ArrayList<>();
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
        plotData.setData(allData);
        if (resProps.getExpMode().equals("cpmg")) {
            plotData.setNames("CPMG", "\u03BD(cpmg)", "R2(\u03BD)");
        } else if (resProps.getExpMode().equals("exp")) {
            plotData.setNames("", "Time (s)", "Intensity");
        } else if (resProps.getExpMode().equals("cest")) {
            plotData.setNames("", "Offset (Hz)", "I(t)/I(0)");
        }

        plotData.setEquations(equations);
        plotData.layoutPlotChildren();
        updateTable(resDatas);
        if (residues != null) {
            updateTableWithPars(mapName, residues, equationName, state, allStates);
            updateEquation(mapName, residues, equationName);
        }

    }

}

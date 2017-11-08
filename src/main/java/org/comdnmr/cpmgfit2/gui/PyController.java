package org.comdnmr.cpmgfit2.gui;

import java.io.File;
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
import javafx.event.ActionEvent;
import javafx.geometry.Side;
import javafx.scene.Node;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.XYChart;
import javafx.scene.control.CheckBox;
import javafx.scene.control.ChoiceBox;
import javafx.scene.control.ContextMenu;
import javafx.scene.control.Label;
import javafx.scene.control.Menu;
import javafx.scene.control.MenuItem;
import javafx.scene.control.Slider;
import javafx.scene.control.TabPane;
import javafx.scene.control.TableColumn;
import javafx.scene.control.TableView;
import javafx.scene.control.cell.PropertyValueFactory;
import javafx.scene.layout.Priority;
import javafx.scene.layout.VBox;
import javafx.scene.paint.Color;
import javafx.scene.shape.Circle;
import javafx.stage.FileChooser;
import javax.script.ScriptException;
import org.comdnmr.cpmgfit2.calc.CPMGFit;
import org.comdnmr.cpmgfit2.calc.CPMGFitResult;
import org.comdnmr.cpmgfit2.calc.CalcRDisp;
import org.comdnmr.cpmgfit2.calc.CurveFit;
import org.comdnmr.cpmgfit2.calc.DataIO;
import org.comdnmr.cpmgfit2.calc.ResidueData;
import org.controlsfx.control.PropertySheet;
import org.comdnmr.cpmgfit2.calc.ParValueInterface;
import org.comdnmr.cpmgfit2.calc.PlotEquation;
import org.comdnmr.cpmgfit2.calc.ProcessingStatus;
import org.comdnmr.cpmgfit2.calc.ResidueFitter;
import org.comdnmr.cpmgfit2.calc.ResidueInfo;
import org.comdnmr.cpmgfit2.calc.ResidueProperties;
import org.controlsfx.control.StatusBar;

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
    TabPane parTabPane;
    @FXML
    PropertySheet propertySheet;

    @FXML
    VBox chartBox;

    @FXML
    ChoiceBox equationSelector;

    @FXML
    Slider kExSlider;

    @FXML
    Slider r2Slider;

    @FXML
    Slider rExSlider;

    @FXML
    Slider pASlider;

    @FXML
    Slider dWSlider;
    @FXML
    Slider field2Slider;

    @FXML
    Label kExValue;

    @FXML
    Label rExValue;

    @FXML
    Label r2Value;
    @FXML
    Label pAValue;
    @FXML
    Label dWValue;
    @FXML
    Label field2Value;
    @FXML
    StatusBar statusBar;
    Circle statusCircle;
    @FXML
    CheckBox absValueModeCheckBox;
    @FXML
    CheckBox nonParBootStrapCheckBox;

    boolean updatingTable = false;

    final ContextMenu axisMenu = new ContextMenu();
    final static double defaultField = 500.0;

    ResidueInfo currentResInfo = null;
    ResidueProperties currentResProps = null;
    ResidueFitter residueFitter;

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
        simSliderAction("");
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
        r2Slider.valueProperty().addListener(e -> {
            simSliderAction("R2");
        });
        kExSlider.valueProperty().addListener(e -> {
            simSliderAction("Kex");
        });
        rExSlider.valueProperty().addListener(e -> {
            simSliderAction("Rex");
        });
        pASlider.valueProperty().addListener(e -> {
            simSliderAction("pA");
        });
        dWSlider.valueProperty().addListener(e -> {
            simSliderAction("dW");
        });
        field2Slider.valueProperty().addListener(e -> {
            simSliderAction("Field2");
        });
        equationSelector.valueProperty().addListener(e -> {
            equationAction();
        });
        residueFitter = new ResidueFitter(this::updateFitProgress, this::updateStatus);
        statusCircle = new Circle(10);
        statusBar.getLeftItems().add(statusCircle);

    }

    void equationAction() {
        String equationName = equationSelector.getValue().toString();
        switch (equationName) {
            case "NOEX":
                r2Slider.setDisable(false);
                rExSlider.setDisable(true);
                kExSlider.setDisable(true);
                pASlider.setDisable(true);
                dWSlider.setDisable(true);
                break;
            case "CPMGFAST":
                r2Slider.setDisable(false);
                rExSlider.setDisable(false);
                kExSlider.setDisable(false);
                pASlider.setDisable(true);
                dWSlider.setDisable(true);
                break;
            case "CPMGSLOW":
                r2Slider.setDisable(false);
                rExSlider.setDisable(true);
                kExSlider.setDisable(false);
                pASlider.setDisable(false);
                dWSlider.setDisable(false);
                break;
            default:
                return;
        }
        // simSliderAction(equationName);

    }

    void simSliderAction(String label) {
        if (updatingTable) {
            return;
        }
        String equationName = equationSelector.getValue().toString();
        if (equationName.equals("CPMGSLOW") && label.equals("Rex")) {
            return;
        }
        double r2 = r2Slider.getValue();
        double kEx = kExSlider.getValue();
        double rEx = rExSlider.getValue();
        double pA = pASlider.getValue();
        double dW = dWSlider.getValue();
        double field2 = field2Slider.getValue();
        r2Value.setText(String.format("%.1f", r2));
        kExValue.setText(String.format("%.1f", kEx));
        rExValue.setText(String.format("%.1f", rEx));
        pAValue.setText(String.format("%.3f", pA));
        dWValue.setText(String.format("%.1f", dW));
        field2Value.setText(String.format("%.0f", field2));
        double[] pars;
        switch (equationName) {
            case "NOEX":
                pars = new double[1];
                pars[0] = r2;

                break;
            case "CPMGFAST":
                pars = new double[3];
                pars[0] = kEx;
                pars[1] = r2;
                pars[2] = rEx;
                break;
            case "CPMGSLOW":
                pars = new double[4];
                pars[0] = kEx;
                pars[1] = pA;
                pars[2] = r2;
                pars[3] = dW;
                int[] map = {0};
                //rEx = CalcRDisp.CPMGEquation.CPMGSLOW.getRex(pars, map);
                rEx = 0.0;
                rExValue.setText(String.format("%.1f", rEx));
                rExSlider.setValue(rEx);
                break;
            default:
                return;
        }
        double[] errs = new double[pars.length];
        int nFields = field2 > (defaultField + 10) ? 2 : 1; // add 10.0 to make sure slider set near to bottom gives 1 field
        double[] fields = new double[nFields];
        fields[0] = 1.0;
        if (nFields > 1) {
            fields[1] = field2 / defaultField;
        }
        updateChartEquations(equationName, pars, errs, fields);
    }

    public void updateChartEquations(String equationName, double[] pars, double[] errs, double[] fields) {
        List<PlotEquation> equations = new ArrayList<>();
        for (int i = 0; i < fields.length; i++) {
            double[] extras = {fields[i] / fields[0]};
            PlotEquation plotEquation = new PlotEquation(equationName, pars, errs, extras);
            equations.add(plotEquation);
        }
        cpmgchart.setEquations(equations);
        cpmgchart.layoutPlotChildren();

    }

    public XYBarChart getActiveChart() {
        return activeChart;
    }

    public XYBarChart addChart() {
        XYBarChart newChart = new XYBarChart();
        int nChildren = chartBox.getChildren().size();
        chartBox.getChildren().add(1, newChart);
        VBox.setVgrow(newChart, Priority.ALWAYS);
        activeChart = newChart;
        return newChart;

    }

    public void updateTable(List<ResidueData> resDatas) {
        ObservableList<ResidueData.DataValue> data = FXCollections.observableArrayList();
        for (ResidueData resData : resDatas) {
            data.addAll(resData.getDataValues());
        }
        resInfoTable.itemsProperty().setValue(data);

        TableColumn<ResidueData.DataValue, String> nameColumn = new TableColumn<>("Name");
        TableColumn<ResidueData.DataValue, String> resColumn = new TableColumn<>("Residue");
        TableColumn<ResidueData.DataValue, Double> xColumn = new TableColumn<>("Vcpmg");
        TableColumn<ResidueData.DataValue, Double> yColumn = new TableColumn<>("Reff");
        TableColumn<ResidueData.DataValue, String> errColumn = new TableColumn<>("Error");
        TableColumn<ResidueData.DataValue, String> peakColumn = new TableColumn<>("Peak");

        nameColumn.setCellValueFactory(new PropertyValueFactory<>("Name"));
        resColumn.setCellValueFactory(new PropertyValueFactory<>("Residue"));
        xColumn.setCellValueFactory(new PropertyValueFactory<>("X"));
        yColumn.setCellValueFactory(new PropertyValueFactory<>("Y"));
        errColumn.setCellValueFactory(new PropertyValueFactory<>("Error"));
        peakColumn.setCellValueFactory(new PropertyValueFactory<>("Peak"));

        resInfoTable.getColumns().clear();
        resInfoTable.getColumns().addAll(nameColumn, resColumn, xColumn, yColumn, errColumn, peakColumn);
    }

    public void updateTableWithPars(String mapName, String[] residues, String equationName, String state) {
        List<ParValueInterface> allParValues = new ArrayList<>();
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
                    updateSliders(parValues, useEquationName);
                }

                allParValues.addAll(parValues);
            }
        }
        updateTableWithPars(allParValues);

    }

    public void updateTableWithPars(ResidueInfo resInfo, String equationName, String state) {
        System.out.println("update table with pars " + equationName + " " + state);
        final String useEquationName;
        if (equationName.equals("best")) {
            useEquationName = resInfo.getBestEquationName();
        } else {
            useEquationName = equationName;
        }
        List<ParValueInterface> parValues = resInfo.getParValues(useEquationName, state);
        currentResInfo = resInfo;
        updateTableWithPars(parValues);

        updateSliders(parValues, useEquationName);
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

    public void updateSliders(List<ParValueInterface> parValues, String equationName) {
        updatingTable = true;
        for (ParValueInterface parValue : parValues) {
            String parName = parValue.getName();
            switch (parName) {
                case "R2":
                    double r2 = parValue.getValue();
                    r2Slider.setValue(r2);
                    r2Value.setText(String.format("%.1f", r2));
                    break;
                case "Rex":
                    double rEx = parValue.getValue();
                    rExSlider.setValue(rEx);
                    rExValue.setText(String.format("%.1f", rEx));
                    break;
                case "Kex":
                    double kEx = parValue.getValue();
                    kExSlider.setValue(kEx);
                    kExValue.setText(String.format("%.1f", kEx));
                    break;
                case "dW":
                    double dW = parValue.getValue();
                    dWSlider.setValue(dW);
                    dWValue.setText(String.format("%.1f", dW));
                    break;
                case "pA":
                    double pA = parValue.getValue();
                    pASlider.setValue(pA);
                    pAValue.setText(String.format("%.3f", pA));
                    break;
                default:
            }
        }
        equationSelector.setValue(equationName);
        updatingTable = false;
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

    public void clearChart() {
        activeChart.getData().clear();
    }

    public void setYAxisType(String setName, String eqnName, String state, String typeName) {
        ObservableList<XYChart.Series<Double, Double>> data = ChartUtil.getParMapData(setName, eqnName, state, typeName);
        XYBarChart chart = activeChart;
        ((XYChart) chart).getData().addAll(data);
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
    }

    void makeAxisMenu() {
        axisMenu.getItems().clear();
        //Residue	 Peak	GrpSz	Group	Equation	   RMS	   AIC	Best	     R2	  R2.sd	    Rex	 Rex.sd	    Kex	 Kex.sd	     pA	  pA.sd	     dW	  dW.sd

        String[] parTypes = {"R2", "Rex", "Kex", "pA", "dW", "RMS", "AIC", "Equation"};
        Map<String, ResidueProperties> residueProps = ChartUtil.residueProperties;
        MenuItem clearItem = new MenuItem("Clear");
        clearItem.setOnAction(e -> clearChart());
        axisMenu.getItems().add(clearItem);
        for (String parType : parTypes) {
            Menu cascade = new Menu(parType);
            axisMenu.getItems().add(cascade);
            System.out.println(parType);
            for (ResidueProperties residueProp : residueProps.values()) {
                String setName = residueProp.getName();
                System.out.println(setName);
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
                            String[] validPars = CalcRDisp.CPMGEquation.valueOf(equationName).getParNames();
                            for (String validPar : validPars) {
                                if (validPar.equals(parType)) {
                                    isValidPar = true;
                                    break;
                                }
                            }
                        }
                        if (isValidPar) {
                            System.out.println("eq " + equationName);
                            Menu cascade3 = new Menu(equationName);
                            cascade2.getItems().add(cascade3);
                            residueProp.getStateStrings().stream().forEach(state -> {
                                System.out.println("ss " + state);
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
    public void fitEquation(ActionEvent event) {
        CPMGFit cpmgFit = new CPMGFit();
        String[] resNums = {String.valueOf(currentResInfo.getResNum())};
        cpmgFit.setData(currentResProps, resNums);
        String equationName = (String) equationSelector.getValue();
        CPMGFitResult fitResult = cpmgFit.doFit(currentResProps, (String) equationSelector.getValue(), absValueModeCheckBox.isSelected(), nonParBootStrapCheckBox.isSelected());
        updateAfterFit(fitResult);
    }

    public void updateAfterFit(CPMGFitResult fitResult) {
        CurveFit curveFit = fitResult.getCurveFit(0);

        List<ParValueInterface> parValues = new ArrayList<>();
        updateTableWithPars(parValues);
        updateSliders(parValues, fitResult.getEquationName());
        double[] pars = curveFit.getEquation().getPars();
        double[] errs = curveFit.getEquation().getErrs();
        updateChartEquations(fitResult.getEquationName(), pars, errs, curveFit.getEquation().getExtras());

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
            residueFitter.fitResidues(currentResProps, allResidues);
        }
    }
    
    public void refreshFit() {
        XYBarChart chart = getActiveChart();
        chart.showInfo();
    }

    @FXML
    public void haltFit(ActionEvent event) {
        residueFitter.haltFit();
    }

    @FXML
    public void saveParameters(ActionEvent event) {
        FileChooser fileChooser = new FileChooser();
        File file = fileChooser.showSaveDialog(axisMenu);
        if (file != null) {
            DataIO.saveParametersToFile(file.getAbsolutePath(), currentResProps);
        }
    }

    public Double updateFitProgress(Double f) {
        if (Platform.isFxApplicationThread()) {
            clearChart();
            statusBar.setProgress(f);
            setYAxisType(currentResProps.getName(), "best", "0:0:0", "RMS");
        } else {
            Platform.runLater(() -> {
                clearChart();
                setYAxisType(currentResProps.getName(), "best", "0:0:0", "RMS");
                statusBar.setProgress(f);
            });
        }
        return null;

    }
    
    public void processingDone() {
        
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

}

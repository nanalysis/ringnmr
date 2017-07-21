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
import javafx.scene.chart.Chart;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.XYChart;
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
                rEx = CalcRDisp.CPMGEquation.CPMGSLOW.getRex(pars);
                rExValue.setText(String.format("%.1f", rEx));
                rExSlider.setValue(rEx);
                break;
            default:
                return;
        }
        int nFields = field2 > (defaultField + 10) ? 2 : 1; // add 10.0 to make sure slider set near to bottom gives 1 field
        double[] fields = new double[nFields];
        fields[0] = 1.0;
        if (nFields > 1) {
            fields[1] = field2 / defaultField;
        }
        updateChartEquations(equationName, pars, fields);
    }

    public void updateChartEquations(String equationName, double[] pars, double[] fields) {
        List<PlotEquation> equations = new ArrayList<>();
        for (int i = 0; i < fields.length; i++) {
            double[] extras = {fields[i] / fields[0]};
            PlotEquation plotEquation = new PlotEquation(equationName, pars, extras);
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

    public void updateTable(ResidueData resData) {
        ObservableList<ResidueData.DataValue> data = FXCollections.observableArrayList();
        data.addAll(resData.getDataValues());
        resInfoTable.itemsProperty().setValue(data);

        TableColumn<ResidueData.DataValue, Double> xColumn = new TableColumn<>("Vcpmg");
        TableColumn<ResidueData.DataValue, Double> yColumn = new TableColumn<>("Reff");
        TableColumn<ResidueData.DataValue, String> peakColumn = new TableColumn<>("Peak");

        xColumn.setCellValueFactory(new PropertyValueFactory<>("X"));
        yColumn.setCellValueFactory(new PropertyValueFactory<>("Y"));
        peakColumn.setCellValueFactory(new PropertyValueFactory<>("Peak"));

        resInfoTable.getColumns().clear();
        resInfoTable.getColumns().addAll(xColumn, yColumn, peakColumn);
    }

    public void updateTableWithPars(ResidueInfo resInfo, String equationName) {
        List<ParValueInterface> parValues = resInfo.getParValues(equationName);
        currentResInfo = resInfo;
        updateTableWithPars(parValues);
        final String useEquationName;
        if (equationName.equals("best")) {
            useEquationName = resInfo.getBestEquationName();
        } else {
            useEquationName = equationName;
        }

        updateSliders(parValues, useEquationName);
    }

    public void updateTableWithPars(List<ParValueInterface> parValues) {
        ObservableList<ParValueInterface> data = FXCollections.observableArrayList();
        data.addAll(parValues);
        parameterTable.itemsProperty().setValue(data);

        TableColumn<ParValueInterface, String> nameColumn = new TableColumn<>("Name");
        TableColumn<ParValueInterface, Double> valueColumn = new TableColumn<>("Value");
        TableColumn<ParValueInterface, Double> errorColumn = new TableColumn<>("Error");

        nameColumn.setCellValueFactory(new PropertyValueFactory<>("Name"));
        valueColumn.setCellValueFactory(new PropertyValueFactory<>("Value"));
        errorColumn.setCellValueFactory(new PropertyValueFactory<>("Error"));

        parameterTable.getColumns().clear();
        parameterTable.getColumns().addAll(nameColumn, valueColumn, errorColumn);
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

    public void selectTableRow(int iRow) {
        parTabPane.getSelectionModel().select(1);
        resInfoTable.getSelectionModel().clearAndSelect(iRow);
        resInfoTable.scrollTo(iRow);
    }

    public void clearChart() {
        activeChart.getData().clear();
    }

    public void setYAxisType(String setName, String eqnName, String typeName) {
        ObservableList<XYChart.Series<Double, Double>> data = ChartUtil.getParMapData(setName, eqnName, typeName);
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
            for (ResidueProperties residueProp : residueProps.values()) {
                String setName = residueProp.getName();
                Menu cascade2 = new Menu(setName);
                cascade.getItems().add(cascade2);
                ArrayList<String> equationNames = new ArrayList<>();
                equationNames.add("best");
                equationNames.addAll(residueProp.getEquationNames());
                for (String equationName : equationNames) {
                    MenuItem cmItem1 = new MenuItem(equationName);
                    cmItem1.setOnAction(e -> setYAxisType(setName, equationName, parType));
                    cascade2.getItems().add(cmItem1);
                }
            }
        }
    }

    @FXML
    public void fitEquation(ActionEvent event) {
        CPMGFit cpmgFit = new CPMGFit();
        String[] resNums = {String.valueOf(currentResInfo.getResNum())};
        cpmgFit.setData(currentResProps.getExperimentData(), resNums);
        String equationName = (String) equationSelector.getValue();
        CPMGFitResult fitResult = cpmgFit.doFit((String) equationSelector.getValue());
        updateAfterFit(fitResult);
    }

    public void updateAfterFit(CPMGFitResult fitResult) {
        List<ParValueInterface> parValues = fitResult.getParValues(0);
        double[] pars = new double[parValues.size()];
        for (int i = 0; i < pars.length; i++) {
            pars[i] = parValues.get(i).getValue();
        }
        updateTableWithPars(parValues);
        updateSliders(parValues, fitResult.getEquationName());
        updateChartEquations(fitResult.getEquationName(), pars, fitResult.getUsedFields());

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
            setYAxisType(currentResProps.getName(), "best", "RMS");
        } else {
            Platform.runLater(() -> {
                clearChart();
                setYAxisType(currentResProps.getName(), "best", "RMS");
                statusBar.setProgress(f);
            });
        }
        return null;

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

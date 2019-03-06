package org.comdnmr.fit.gui;

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
import javafx.embed.swing.SwingFXUtils;
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
import java.util.Optional;
import javafx.beans.property.SimpleStringProperty;
import java.util.Random;
import javafx.print.PrinterJob;
import javafx.scene.canvas.Canvas;
import javafx.scene.canvas.GraphicsContext;
import javafx.scene.control.Alert;
import javafx.scene.control.SplitPane;
import javafx.scene.control.Label;
import javafx.scene.control.TextField;
import javafx.scene.control.Button;
import javafx.scene.control.ContentDisplay;
import javafx.scene.control.ToolBar;
import javafx.scene.input.MouseEvent;
import javafx.scene.layout.Pane;
import org.comdnmr.fit.calc.ExperimentData.Nuclei;
import org.comdnmr.fit.calc.FitModel;
import org.comdnmr.fit.calc.R1RhoFit;
import static org.comdnmr.fit.gui.MainApp.preferencesController;
import static org.comdnmr.fit.gui.MainApp.console;
import static org.comdnmr.fit.gui.MainApp.primaryStage;
import org.comdnmr.utils.NMRFxClient;
import org.controlsfx.dialog.ExceptionDialog;
import org.nmrfx.chart.Axis;
import org.nmrfx.chart.DataSeries;
import org.nmrfx.chart.XYEValue;
import org.nmrfx.chart.XYValue;
import org.nmrfx.graphicsio.GraphicsIOException;
import org.nmrfx.graphicsio.SVGGraphicsContext;

public class PyController implements Initializable {

    public static PyController mainController;
    ResidueChart activeChart;
    List<ResidueChart> barCharts = new ArrayList<>();

    SSPainter ssPainter = null;
    @FXML
    SSRegion ssregion;

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
    NMRFxClient cl;

    //final ContextMenu axisMenu = new ContextMenu();
    static double defaultField = 500.0;

    ResidueInfo currentResInfo = null;
    ResidueProperties currentResProps = null;
    ResidueFitter residueFitter;
    String currentMapName = "";
    String[] currentResidues;
    List<String> fittingResidues = new ArrayList<>();
    String currentState = "";
    String currentEquationName;
    List<int[]> currentStates = new ArrayList<>();

    boolean simulate = true;

    CPMGFitResult fitResult;
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
        } else if (getFittingMode().equals("cest")) {
            simControls = new CESTControls();
            xLowerBoundTextField.setText("-20.0");
            xUpperBoundTextField.setText("20.0");
            if (currentResProps.getExperimentData() != null) {
                double[] xVals = currentResProps.getExperimentData().stream().findFirst().get().getXVals();
                xLowerBoundTextField.setText(String.valueOf(Math.floor(xVals[1] / 2) * 2));
                xUpperBoundTextField.setText(String.valueOf(Math.ceil(xVals[xVals.length - 1] / 2) * 2));
            }
            yLowerBoundTextField.setText("0.0");
            yUpperBoundTextField.setText("1.0");
            xTickTextField.setText("1.0");
            yTickTextField.setText("0.25");
        } else if (getFittingMode().equals("exp")) {
            simControls = new ExpControls();
            xLowerBoundTextField.setText("0.0");
            xUpperBoundTextField.setText("1.25");
            yLowerBoundTextField.setText("0.0");
            yUpperBoundTextField.setText("100.0");
            xTickTextField.setText("0.25");
            yTickTextField.setText("10.0");
        } else if (getFittingMode().equals("r1rho")) {
            simControls = new R1RhoControls();
            xLowerBoundTextField.setText("-20.0");
            xUpperBoundTextField.setText("20.0");
            if (currentResProps.getExperimentData() != null) {
                double[] xVals = currentResProps.getExperimentData().stream().findFirst().get().getXVals();
                xLowerBoundTextField.setText(String.valueOf(Math.floor(xVals[1] / 2) * 2));
                xUpperBoundTextField.setText(String.valueOf(Math.ceil(xVals[xVals.length - 1] / 2) * 2));
            }
            yLowerBoundTextField.setText("0.0");
            yUpperBoundTextField.setText("50.0");
            xTickTextField.setText("1.0");
            yTickTextField.setText("0.25");
        } else {
            System.out.println("Error: no fitting mode selected.");
        }
        VBox vBox = simControls.makeControls(mainController);
        simPane.centerProperty().set(vBox);
        residueFitter = new ResidueFitter(this::updateFitProgress, this::updateStatus);
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
        calcErrorsCheckBox.selectedProperty().addListener(e -> FitModel.setCalcError(calcErrorsCheckBox.isSelected()));
        calcErrorsCheckBox.setSelected(true);
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
            } else {
                if ((x > xAxis.getXOrigin()) && (x < xAxis.getXOrigin() + xAxis.getWidth())) {
                    if ((y < yAxis.getYOrigin()) && (y > xAxis.getYOrigin() - yAxis.getHeight())) {
                        activeChart = residueChart;
                        break;
                    }
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
        updateXYChartLabels();
    }

    public void setSimControls() {
        boolean update = false;
        if (getSimMode().equals("cpmg") && !(simControls instanceof CPMGControls)) {
            simControls = new CPMGControls();
            update = true;
        } else if (getSimMode().equals("exp") && !(simControls instanceof ExpControls)) {
            simControls = new ExpControls();
            update = true;
        } else if (getSimMode().equals("cest") && !(simControls instanceof CESTControls)) {
            simControls = new CESTControls();
            ((CESTControls) simControls).updateDeltaLimits();
            update = true;
        } else if (getFittingMode().equals("r1rho") && !(simControls instanceof R1RhoControls)) {
            simControls = new R1RhoControls();
            ((R1RhoControls) simControls).updateDeltaLimits();
            update = true;
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
            simControls.updateEquations(equationChoice, CPMGFit.getEquationNames());
        } else if (mode.equals("exp")) {
            simControls.updateEquations(equationChoice, ExpFit.getEquationNames());
        } else if (mode.equals("cest")) {
            simControls.updateEquations(equationChoice, CESTFit.getEquationNames());
        } else if (mode.equals("r1rho")) {
            simControls.updateEquations(equationChoice, R1RhoFit.getEquationNames());
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
            ResidueChart chart = getActiveChart();
            chart.showInfo(chart.currentSeriesName, 0, res, false);
        }
    }

    public void firstResidue(ActionEvent event) {
        if (currentResidues != null) {
            int res = ChartUtil.minRes;
            ResidueChart chart = getActiveChart();
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
            ResidueChart chart = getActiveChart();
            chart.showInfo(chart.currentSeriesName, 0, res, false);
        }
    }

    public void lastResidue(ActionEvent event) {
        if (currentResidues != null) {
            int res = ChartUtil.maxRes;
            ResidueChart chart = getActiveChart();
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
    public void loadSecondaryStructure() {
        List<SecondaryStructure> ssValues = SecondaryStructure.loadFromFile();
        if (!ssValues.isEmpty()) {
            ssPainter = new SSPainter(barPlotCanvas, ssValues);
        }
        resizeBarPlotCanvas();
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
            xychart.setBounds(0.0, 1.25, 0.0, 1.25, 0.25, 0.25);
            xLowerBoundTextField.setText("0.0");
            xUpperBoundTextField.setText("1.25");
            yLowerBoundTextField.setText("0.0");
            yUpperBoundTextField.setText("100.0");
            xTickTextField.setText("0.25");
            yTickTextField.setText("10.0");
        } else if ((simControls instanceof CESTControls)) {
            xychart.setNames("CEST", "Offset (PPM)", "I(t)/I(0)", "20");
            xychart.setBounds(-20, 20, 0.0, 1.0, 2.0, 0.25);
            xLowerBoundTextField.setText("-20.0");
            xUpperBoundTextField.setText("20.0");
            if (currentResProps != null) {
                double[] xVals = currentResProps.getExperimentData().stream().findFirst().get().getXVals();
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
            xychart.setBounds(-20, 20, 0.0, 1.0, 2.0, 0.25);
            xLowerBoundTextField.setText("-20.0");
            xUpperBoundTextField.setText("20.0");
            if (currentResProps != null) {
                double[] xVals = currentResProps.getExperimentData().stream().findFirst().get().getXVals();
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
        xLowerBoundTextField.setText(Double.toString(bounds[0]));
        xUpperBoundTextField.setText(Double.toString(bounds[1]));
        yLowerBoundTextField.setText(Double.toString(bounds[2]));
        yUpperBoundTextField.setText(Double.toString(bounds[3]));
        if (simControls instanceof CESTControls) {
            ((CESTControls) simControls).updateDeltaLimits(bounds[0], bounds[1]);
        }

    }

    public void updateChartEquations(String equationName, double[] pars, double[] errs, double[] fields) {
        List<GUIPlotEquation> equations = new ArrayList<>();
        for (int i = 0; i < fields.length; i++) {
            double[] extras = {fields[i] / fields[0]};
            //System.out.println("updateChartEquations got called with extras length = "+extras.length);
            GUIPlotEquation plotEquation = new GUIPlotEquation(equationName, pars, errs, extras);
            equations.add(plotEquation);
        }
        showEquations(equations);
    }

    public void showEquations(List<GUIPlotEquation> equations) {
        xychart.setEquations(equations);
//        Optional<Double> rms = rms();

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

    public void updateTable(List<ResidueData> resDatas) {
        ObservableList<ResidueData.DataValue> data = FXCollections.observableArrayList();
        for (ResidueData resData : resDatas) {
            data.addAll(resData.getDataValues());
        }
        resInfoTable.itemsProperty().setValue(data);

        TableColumn<ResidueData.DataValue, String> nameColumn = new TableColumn<>("Name");
        TableColumn<ResidueData.DataValue, String> resColumn = new TableColumn<>("Residue");
        TableColumn<ResidueData.DataValue, String> errColumn = new TableColumn<>("Error");
        TableColumn<ResidueData.DataValue, String> peakColumn = new TableColumn<>("Peak");

        nameColumn.setCellValueFactory(new PropertyValueFactory<>("Name"));
        resColumn.setCellValueFactory(new PropertyValueFactory<>("Residue"));
        errColumn.setCellValueFactory(new PropertyValueFactory<>("Error"));
        peakColumn.setCellValueFactory(new PropertyValueFactory<>("Peak"));

        resInfoTable.getColumns().clear();
        resInfoTable.getColumns().addAll(nameColumn, resColumn, errColumn, peakColumn);

        if (getFittingMode().equals("cpmg")) {
            TableColumn<ResidueData.DataValue, Double> xColumn = new TableColumn<>("Vcpmg");
            TableColumn<ResidueData.DataValue, Double> yColumn = new TableColumn<>("Reff");

            xColumn.setCellValueFactory(new PropertyValueFactory<>("X0"));
            yColumn.setCellValueFactory(new PropertyValueFactory<>("Y"));

            resInfoTable.getColumns().clear();
            resInfoTable.getColumns().addAll(nameColumn, resColumn, xColumn, yColumn, errColumn, peakColumn);
        } else if (getFittingMode().equals("cest")) {
            TableColumn<ResidueData.DataValue, Double> x0Column = new TableColumn<>("Offset");
            TableColumn<ResidueData.DataValue, Double> x1Column = new TableColumn<>("B1 Field");
            TableColumn<ResidueData.DataValue, Double> yColumn = new TableColumn<>("Intensity");

            x0Column.setCellValueFactory(new PropertyValueFactory<>("X0"));
            x1Column.setCellValueFactory(new PropertyValueFactory<>("X1"));
            yColumn.setCellValueFactory(new PropertyValueFactory<>("Y"));

            resInfoTable.getColumns().clear();
            resInfoTable.getColumns().addAll(nameColumn, resColumn, x0Column, x1Column, yColumn, errColumn, peakColumn);
        } else if (getFittingMode().equals("r1rho")) {
            TableColumn<ResidueData.DataValue, Double> x0Column = new TableColumn<>("Offset");
            TableColumn<ResidueData.DataValue, Double> x1Column = new TableColumn<>("B1 Field");
            TableColumn<ResidueData.DataValue, Double> yColumn = new TableColumn<>("Intensity");

            x0Column.setCellValueFactory(new PropertyValueFactory<>("X0"));
            x1Column.setCellValueFactory(new PropertyValueFactory<>("X1"));
            yColumn.setCellValueFactory(new PropertyValueFactory<>("Y"));

            resInfoTable.getColumns().clear();
            resInfoTable.getColumns().addAll(nameColumn, resColumn, x0Column, x1Column, yColumn, errColumn, peakColumn);
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
                        String rChiSq = String.format("%.2f", curveSet.getParMap().get("rChiSq"));
                        aicLabel.setText(aic);
                        rmsLabel.setText(rms);
                        rChiSqLabel.setText(rChiSq);
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
        ObservableList<ParValueInterface> data = FXCollections.observableArrayList();

        for (ParValueInterface parValue : parValues) {
            String state = parValue.getState();
            String parName = parValue.getName();
            boolean ok = true;
            // remove redundant parameters 
            if ((state != null) && (parName != null)) {
                if (!state.equals("0:0:0") && (parName.equals("Kex") || parName.equals("dPPM")
                        || parName.equals("pA") || parName.equals("Rex"))) {
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

    public String getPeakNumFromTable() { //getPeakNumFromTable(String seriesName, int index)
        parTabPane.getSelectionModel().select(1);
        List<ResidueData.DataValue> data = resInfoTable.getItems();
        int iRow = 0;
        int peakNum = 0;
        String name = "";
//        for (ResidueData.DataValue dValue : data) {
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

    public void setYAxisType(String setName, String eqnName, String state, String typeName) {
        ObservableList<DataSeries> data = ChartUtil.getParMapData(setName, eqnName, state, typeName);
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

        if (typeName.equals("Kex")) {
            chart.yAxis.setLabel("Kex");
//            chart.setUnifyYAxes(true); // fixme
        } else {
            chart.yAxis.setLabel(typeName);
//            chart.setUnifyYAxes(false); // fixme
        }
        chart.setTitle(setName);
        if (chart.getData().size() > 1) {
//            chart.setLegendVisible(true); // fixme
//            chart.setLegendSide(Side.TOP);
        } else {
            // chart.setLegendVisible(false); //fixme
        }
        currentResProps = ChartUtil.residueProperties.get(setName);
        chart.setResProps(currentResProps);
        refreshResidueCharts();
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
                List<ParValueInterface> guesses = equationFitter.guessPars(equationName);
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
                fitResult = equationFitter.doFit(equationName, sliderGuesses);
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
            //System.out.println("Fit button residueProperties = " + residueProperties);
            //System.out.println("Fit button expData = " + residueProps.getExperimentData("cest"));
            Optional<ExperimentData> optionalData = Optional.empty();
            if (currentResProps != null) {
                optionalData = currentResProps.getExperimentData().stream().findFirst();
            }

            if (optionalData.isPresent() && optionalData.get().getExtras().size() > 0) {
                for (ExperimentData expData : currentResProps.getExperimentData()) {
                    double[] pars = curveFit.getEquation().getPars(); //pars = getPars(equationName);
                    double[] errs = curveFit.getEquation().getErrs(); //double[] errs = new double[pars.length];
                    double[] extras = new double[3];
                    extras[0] = expData.getNucleusField();
                    extras[1] = expData.getExtras().get(0);
                    extras[2] = expData.getExtras().get(1);
//                    System.out.println("Fit button expData extras size = " + expData.getExtras().size() + " extra[1] = " + extras[1]);
                    GUIPlotEquation plotEquation = new GUIPlotEquation(equationName, pars, errs, extras);
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
//                    extras[0] = expData.getField();
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
                for (int i = 0; i < simExtras.length; i++) {
//                    System.out.println("simextras " + i + " " + simExtras[i]);
//                    extras[i + 1] = simExtras[i];
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
                GUIPlotEquation plotEquation = new GUIPlotEquation(equationName, pars, errs, extras);

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
        fittingResidues.clear();
        fittingResidues.addAll(ResidueChart.selectedResidues);
        groupResidues.addAll(ResidueChart.selectedResidues);
        if (!groupResidues.isEmpty()) {
            allResidues.add(groupResidues);
            residueFitter.fitResidues(currentResProps, allResidues);
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
        File file = fileChooser.showSaveDialog(MainApp.primaryStage);
        if (file != null) {
            DataIO.saveResultsFile(file.getAbsolutePath(), currentResProps, false);
        }
    }

    public Double updateFitProgress(Double f) {
        ResidueChart chart = getActiveChart();
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
        } else if (getFittingMode().equals("r1rho")) {
            return new R1RhoFit();
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
        String[] r1rhoTypes = {"kex", "pb", "deltaA0", "deltaB0", "R1A", "R1B", "R2A", "R2B", "RMS", "AIC", "Equation"};
        String[] nullTypes = {"RMS", "AIC", "Equation"};
        if (getFittingMode().equals("exp")) {
            return expTypes;
        } else if (getFittingMode().equals("cpmg")) {
            return cpmgTypes;
        } else if (getFittingMode().equals("cest")) {
            return cestTypes;
        } else if (getFittingMode().equals("r1rho")) {
            return r1rhoTypes;
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
        // fixme xychart.clear();
        barCharts.remove(activeChart);
        activeChart = barCharts.get(0);
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
        ArrayList<GUIPlotEquation> equations = new ArrayList<>();
        ObservableList<DataSeries> allData = FXCollections.observableArrayList();
        List<ResidueData> resDatas = new ArrayList<>();
        List<int[]> allStates = new ArrayList<>();
        if ((resProps != null) && (residues != null)) {
            double maxY = getMaxY(resProps, equationName, mapName, state, residues);
            System.out.println("max Y " + maxY);
            int iSeries = 0;
            for (ExperimentData expData : resProps.getExperimentData()) {
                if (!ResidueProperties.matchStateString(state, expData.getState())) {
                    continue;
                }
                String expName = expData.getName();
                for (String resNum : residues) {
                    resDatas.add(expData.getResidueData(resNum));
                    DataSeries series = ChartUtil.getMapData(mapName, expName, resNum);
                    series.setStroke(PlotData.colors[iSeries % 8]);
                    series.setFill(PlotData.colors[iSeries % 8]);

                    allData.add(series);
                    GUIPlotEquation equation = ChartUtil.getEquation(expData,
                            mapName, resNum, equationName, expData.getState(),
                            expData.getNucleusField());
                    equation.setScaleValue(maxY);
                    equation.setColor(PlotData.colors[iSeries % 8]);
                    if (equation != null) {
                        equations.add(equation);
                    } else {
                        System.out.println("null eq");
                    }
                    iSeries++;
                }

                int[] states = resProps.getStateIndices(0, expData);
                allStates.add(states);
            }
        }
        updateTable(resDatas);
        if (residues != null) {
            updateTableWithPars(mapName, residues, equationName, state, allStates);
            updateEquation(mapName, residues, equationName);
        }
        plotData.setData(allData);
        setBounds();
        plotData.setEquations(equations);
    }

    public double getMaxY(ResidueProperties resProps, String equationName, String mapName, String state, String[] residues) {
        double maxValue = Double.NEGATIVE_INFINITY;
        if ((resProps != null) && (residues != null)) {
            for (ExperimentData expData : resProps.getExperimentData()) {
                if (!ResidueProperties.matchStateString(state, expData.getState())) {
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
                        System.out.println(minX + " " + valueY);
                        if (valueY > maxValue) {
                            maxValue = valueY;
                        }
                    } else {
//                        DataSeries series = ChartUtil.getMapData(mapName, expName, resNum);
//                        double valueY = series.getValues().stream().max(XYValue::getYValue).get();
//                        if (valueY > maxValue) {
//                            maxValue = valueY;
//                        }

                    }
                }
            }
        }
        return maxValue;
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
        System.out.println("xval len " + xValues.size());
        for (int j = 0; j < extras.length; j++) {
            System.out.println("set sim " + j + " " + extras[j]);
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
        List<ParValueInterface> guesses = equationFitter.guessPars(equationName);
        simControls.updateSliders(guesses, equationName);
        simControls.simSliderAction("");
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
        fitResult = equationFitter.doFit(equationName, sliderGuesses);
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

    @FXML
    void showSimData(ActionEvent e) {
        ObservableList<DataSeries> allData = FXCollections.observableArrayList();
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
                double yValue = eqn.calculate(sliderGuesses, ax, eqn.getExtra(0));
                yValue += sdev * rand.nextGaussian();
                XYValue dataPoint = new XYValue(xValue, yValue);
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

}

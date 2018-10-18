package org.comdnmr.fit.gui;

import de.jensd.fx.glyphs.GlyphsDude;
import de.jensd.fx.glyphs.fontawesome.FontAwesomeIcon;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
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
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Optional;
import java.util.Set;
import javafx.beans.property.SimpleStringProperty;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.Random;
import javafx.event.EventHandler;
import javafx.scene.control.Alert;
import javafx.scene.control.TextArea;
import javafx.scene.control.SplitPane;
import static org.comdnmr.fit.gui.ChartUtil.residueProperties;
import javafx.scene.Scene;
import javafx.scene.control.Label;
import javafx.scene.control.TextField;
import javafx.scene.control.Button;
import javafx.scene.control.ContentDisplay;
import javafx.scene.control.ToolBar;
import javafx.scene.input.KeyCode;
import javafx.scene.input.KeyEvent;
import javafx.scene.layout.GridPane;
import static org.comdnmr.fit.calc.DataIO.loadPeakFile;
import static org.comdnmr.fit.calc.DataIO.loadTextFile;
import org.comdnmr.fit.calc.FitModel;
import static org.comdnmr.fit.gui.MainApp.preferencesController;
import static org.comdnmr.fit.gui.MainApp.console;
import static org.comdnmr.fit.gui.MainApp.primaryStage;
import org.controlsfx.control.textfield.TextFields;
import org.python.util.InteractiveInterpreter;
import org.yaml.snakeyaml.Yaml;

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

    GridPane inputInfoDisplay = new GridPane();
    Scene inputScene = new Scene(inputInfoDisplay, 600, 600);
    Stage infoStage = new Stage();
    TextField chosenFileLabel = new TextField();
    TextField chosenXPK2FileLabel = new TextField();
    TextField chosenParamFileLabel = TextFields.createClearableTextField();
    TextField fieldTextField = new TextField();
    TextField tempTextField = new TextField();
    TextField nucTextField = new TextField();
    TextField pTextField = new TextField();
    TextField modeTextField = new TextField();
    TextField tauTextField = new TextField();
    TextArea xValTextArea = new TextArea();
    ChoiceBox<String> fitModeChoice = new ChoiceBox<>();
    TextField B1TextField = new TextField();
    TextField TexTextField = new TextField();
    TextField yamlTextField = new TextField();
    CheckBox ppmBox = new CheckBox("ppm to Hz");
    ChoiceBox<String> errModeChoice = new ChoiceBox<>();
    TextField errPercentTextField = new TextField();
    ArrayList<HashMap<String, Object>> dataList = new ArrayList();
    Button fileChoiceButton = new Button();
    Button xpk2ChoiceButton = new Button();
    Button paramFileChoiceButton = new Button();
    Button addButton = new Button();
    Button clearButton = new Button();
    Button yamlButton = new Button();
    Button loadButton = new Button();

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

    public void inputParameters(ActionEvent event) {

        infoStage.setTitle("Input Data Parameters");
        Label fileLabel = new Label("  File:  ");
        Label xpk2FileLabel = new Label("  XPK2 File:  ");
        Label fitFileLabel = new Label("  Fit Parameter File:  ");
        Label fieldLabel = new Label("  Field:  ");
        Label tempLabel = new Label("  Temperature:  ");
        Label nucLabel = new Label("  Nucleus:  ");
        Label pLabel = new Label("  Pressure:  ");
        Label modeLabel = new Label("  Mode:  ");
        Label tauLabel = new Label("  Tau:  ");
        Label xValLabel = new Label("  X Values:  ");
        Label fitModeLabel = new Label("  Experiment Type:  ");
        Label B1FieldLabel = new Label("  B1 Field:  ");
        Label TexLabel = new Label("  Tex:  ");
        Label yamlLabel = new Label("  YAML File:  ");
        Label errModeLabel = new Label("  Error Mode:  ");
        Label errPercentLabel = new Label("  Error Value:  ");

        Label[] labels = {fitModeLabel, fileLabel, xpk2FileLabel, fitFileLabel, fieldLabel, tempLabel, nucLabel, pLabel, modeLabel, 
            tauLabel, B1FieldLabel, TexLabel, errModeLabel, errPercentLabel, xValLabel, yamlLabel};

//        Button fileChoiceButton = new Button();
        fileChoiceButton.setOnAction(e -> chooseFile(e));
        fileChoiceButton.setText("Browse");
        chosenFileLabel.setText("");

//        Button xpk2ChoiceButton = new Button();
        xpk2ChoiceButton.setOnAction(e -> chooseXPK2File(e));
        xpk2ChoiceButton.setText("Browse");
        chosenXPK2FileLabel.setText("");

//        Button paramFileChoiceButton = new Button();
        paramFileChoiceButton.setOnAction(e -> chooseParamFile(e));
        paramFileChoiceButton.setText("Browse");
        chosenParamFileLabel.setText("");

        ppmBox.setSelected(false);
        
        double textFieldWidth = 100;
        double xValAreaWidth = 150; //240;

        fieldTextField.setText("");
        tempTextField.setText("25.0");
        nucTextField.setText("");
        pTextField.setText("20.0");
        modeTextField.setText("mpk2");
        tauTextField.setText("0.04");
        B1TextField.setText("20.0");
        TexTextField.setText("0.5");
        errPercentTextField.setText("5");
        errPercentTextField.setMaxWidth(textFieldWidth);
        xValTextArea.setText("");
        xValTextArea.setMaxWidth(xValAreaWidth);
        xValTextArea.setWrapText(true);
        yamlTextField.setText("");

        TextField[] texts = {fieldTextField, tempTextField, nucTextField, pTextField, modeTextField, tauTextField, B1TextField, TexTextField};

        inputInfoDisplay.getChildren().clear();

        for (int i = 0; i < labels.length; i++) {
            inputInfoDisplay.add(labels[i], 0, i);
        }
        for (int i = 0; i < texts.length; i++) {
            inputInfoDisplay.add(texts[i], 1, i + 4);
            texts[i].setMaxWidth(textFieldWidth);
        }
        
        fitModeChoice.getItems().add("Select");
        fitModeChoice.getItems().add("CPMG");
        fitModeChoice.getItems().add("EXP");
        fitModeChoice.getItems().add("CEST");
        fitModeChoice.setValue("Select");

        fitModeChoice.valueProperty().addListener(x -> {
            updateInfoInterface();
        });

        EventHandler<ActionEvent> boxevent = new EventHandler<ActionEvent>() {

            public void handle(ActionEvent e) {
                String[] xvals = xValTextArea.getText().split("\t");
                ArrayList<Double> fxvals = new ArrayList();
                String xString = "";
                if (fitModeChoice.getSelectionModel().getSelectedItem().equals("CEST") && ppmBox.isSelected()) {
                    for (int i = 0; i < xvals.length; i++) {
                        fxvals.add(Double.parseDouble(xvals[i]) * Double.parseDouble(fieldTextField.getText()));
                        xString += fxvals.get(i).toString() + "\t";
                    }
                    xValTextArea.setText(xString);
                } else if (fitModeChoice.getSelectionModel().getSelectedItem().equals("CEST") && !ppmBox.isSelected()) {
                    for (int i = 0; i < xvals.length; i++) {
                        fxvals.add(Double.parseDouble(xvals[i]) / Double.parseDouble(fieldTextField.getText()));
                        xString += fxvals.get(i).toString() + "\t";
                    }
                    xValTextArea.setText(xString);
                }
            }

        };

        // set event to checkbox 
        ppmBox.setOnAction(boxevent);
        
        errModeChoice.getItems().add("percent");
        errModeChoice.getItems().add("replicates");
        errModeChoice.getItems().add("noise");
        errModeChoice.setValue("percent");

        inputInfoDisplay.add(fitModeChoice, 1, 0);
        inputInfoDisplay.add(fileChoiceButton, 2, 1);
        inputInfoDisplay.add(chosenFileLabel, 1, 1);
        inputInfoDisplay.add(xpk2ChoiceButton, 2, 2);
        inputInfoDisplay.add(chosenXPK2FileLabel, 1, 2);
        inputInfoDisplay.add(paramFileChoiceButton, 2, 3);
        inputInfoDisplay.add(chosenParamFileLabel, 1, 3);
        inputInfoDisplay.add(errModeChoice, 1, labels.length - 4);
        inputInfoDisplay.add(errPercentTextField, 1, labels.length - 3);
        inputInfoDisplay.add(xValTextArea, 1, labels.length - 2, 1, 1);
        inputInfoDisplay.add(ppmBox, 2, labels.length - 2);
        inputInfoDisplay.add(yamlTextField, 1, labels.length - 1);
        
        chosenFileLabel.setMaxWidth(200);
        chosenXPK2FileLabel.setMaxWidth(200);
        chosenParamFileLabel.setMaxWidth(200);

//        Button addButton = new Button();
        addButton.setOnAction(e -> addInfo(e));
        addButton.setText("Add to Data List");
        inputInfoDisplay.add(addButton, 1, labels.length);

//        Button clearButton = new Button();
        clearButton.setOnAction(e -> clearDataList(e));
        clearButton.setText("Clear Data List");
        inputInfoDisplay.add(clearButton, 1, labels.length + 1);

//        Button yamlButton = new Button();
        yamlButton.setOnAction(e -> makeYAML(e));
        yamlButton.setText("Create YAML");
        inputInfoDisplay.add(yamlButton, 2, labels.length - 1);

//        Button loadButton = new Button();
        loadButton.setOnAction(e -> loadInfo(e));
        loadButton.setText("Load");
        inputInfoDisplay.add(loadButton, 2, labels.length + 1);

        updateInfoInterface();

        infoStage.setScene(inputScene);
        infoStage.show();

    }

    public void updateInfoInterface() { 
        if (fitModeChoice.getSelectionModel().getSelectedItem().equals("Select")) {
            B1TextField.setDisable(true);
            tauTextField.setDisable(true);
            ppmBox.setDisable(true);
            chosenFileLabel.setDisable(true);
            chosenXPK2FileLabel.setDisable(true);
            chosenParamFileLabel.setDisable(true);
            fieldTextField.setDisable(true);
            tempTextField.setDisable(true);
            nucTextField.setDisable(true);
            pTextField.setDisable(true);
            modeTextField.setDisable(true);
            TexTextField.setDisable(true);
            errPercentTextField.setDisable(true);
            xValTextArea.setDisable(true);
            yamlTextField.setDisable(true);
            errModeChoice.setDisable(true);
            fileChoiceButton.setDisable(true);
            xpk2ChoiceButton.setDisable(true);
            paramFileChoiceButton.setDisable(true);
            addButton.setDisable(true);
            clearButton.setDisable(true);
            yamlButton.setDisable(true);
            loadButton.setDisable(true);
        } else if (!fitModeChoice.getSelectionModel().getSelectedItem().equals("Select")) {
            chosenFileLabel.setDisable(false);
            chosenXPK2FileLabel.setDisable(false);
            chosenParamFileLabel.setDisable(false);
            fieldTextField.setDisable(false);
            tempTextField.setDisable(false);
            nucTextField.setDisable(false);
            pTextField.setDisable(false);
            modeTextField.setDisable(false);
            TexTextField.setDisable(false);
            errPercentTextField.setDisable(false);
            xValTextArea.setDisable(false);
            yamlTextField.setDisable(false);
            errModeChoice.setDisable(false);
            fileChoiceButton.setDisable(false);
            xpk2ChoiceButton.setDisable(false);
            paramFileChoiceButton.setDisable(false);
            addButton.setDisable(false);
            clearButton.setDisable(false);
            yamlButton.setDisable(false);
            loadButton.setDisable(false);
            if (fitModeChoice.getSelectionModel().getSelectedItem().equals("CPMG")) {
                B1TextField.setDisable(true);
                tauTextField.setDisable(false);
                ppmBox.setDisable(true);
            } else if (fitModeChoice.getSelectionModel().getSelectedItem().equals("EXP")) {
                B1TextField.setDisable(true);
                tauTextField.setDisable(true);
                ppmBox.setDisable(true);
            } else if (fitModeChoice.getSelectionModel().getSelectedItem().equals("CEST")) {
                B1TextField.setDisable(false);
                tauTextField.setDisable(true);
                ppmBox.setDisable(false);
            }
        }

    }

    public void chooseFile(ActionEvent event) {
        FileChooser fileChooser = new FileChooser();
        fileChooser.getExtensionFilters().addAll(new FileChooser.ExtensionFilter("mpk2 or txt File", "*.mpk2", "*.txt"));
        File file = fileChooser.showOpenDialog(infoStage);
        String directory = file.getParent();
        String fileName = file.getName();
        chosenFileLabel.setText(directory + "/" + fileName);

        Path path = Paths.get(directory + "/" + fileName);

        try (BufferedReader fileReader = Files.newBufferedReader(path)) {
            while (true) {
                String line = fileReader.readLine();
                if (line == null) {
                    break;
                }
                String sline = line.trim();
                if (sline.length() == 0) {
                    continue;
                }
                if (fileName.endsWith(".mpk2") && !sline.startsWith("id") || fileName.endsWith(".txt") && !sline.startsWith("Residue")) {
                    break;
                }
                String xVals = null;
                if (sline.startsWith("id")) {
                    xVals = line.substring(13);
                } else if (sline.startsWith("Residue")) {
                    if (sline.contains("Hz")) {
                        xVals = line.substring(8).replaceAll(" Hz\t", "\t").replaceAll(" Hz", "");
                    } else if (sline.contains("1/S")) {
                        xVals = line.substring(8).replaceAll(" 1/S\t", "\t").replaceAll(" 1/S", "");
                        StringBuilder xVals1 = new StringBuilder(xVals);
                        xVals1.insert(0, "0.0\t");
                        xVals = xVals1.toString();
                    } else if (sline.contains("S")) {
                        xVals = line.substring(8).replaceAll(" S\t", "\t").replaceAll(" S", "");
                    } else if (sline.contains("On")) {
                        xVals = line.substring(8).replaceAll(" On\t", "\t").replaceAll(" On", "");
                    }
                }
                xValTextArea.setText(xVals);
            }
        } catch (IOException ioe) {
            ioe.printStackTrace();
        }
    }

    public void chooseXPK2File(ActionEvent event) {
        FileChooser fileChooser = new FileChooser();
        fileChooser.getExtensionFilters().addAll(new FileChooser.ExtensionFilter("xpk2 File", "*.xpk2"));
        File file = fileChooser.showOpenDialog(infoStage);
        String directory = file.getParent();
        String fileName = file.getName();
        chosenXPK2FileLabel.setText(directory + "/" + fileName);

        Path path1 = Paths.get(directory + "/" + fileName);

        List<String[]> head = new ArrayList<>();

        try (BufferedReader fileReader = Files.newBufferedReader(path1)) {
            while (true) {
                String line = fileReader.readLine();
                if (line == null) {
                    break;
                }
                String sline = line.trim();
                if (sline.length() == 0) {
                    continue;
                }
                if (sline.startsWith("id")) {
                    break;
                }
                String[] sline1 = line.split("\t", -1);
                head.add(sline1);
            }
        } catch (IOException ioe) {
            ioe.printStackTrace();
        }
        int sfInd = Arrays.asList(head.get(2)).indexOf("sf");
        int codeInd = Arrays.asList(head.get(2)).indexOf("code");
        String field = Arrays.asList(head.get(4)).get(sfInd);
        String nuc = Arrays.asList(head.get(4)).get(codeInd);
        nuc = nuc.replaceAll("[^a-zA-Z]", "");
        nucTextField.setText(nuc);
        fieldTextField.setText(field);
    }

    public void chooseParamFile(ActionEvent event) {
        FileChooser fileChooser = new FileChooser();
        fileChooser.getExtensionFilters().addAll(new FileChooser.ExtensionFilter("txt File", "*.txt"));
        File file = fileChooser.showOpenDialog(infoStage);
        String directory = file.getParent();
        String fileName = file.getName();
        chosenParamFileLabel.setText(directory + "/" + fileName);
    }

    public void resetParamFile(ActionEvent event) {
        chosenParamFileLabel.setText("");
    }

    public void addInfo(ActionEvent event) {
        addInfo();
    }

    public void addInfo() {
        HashMap hm = new HashMap();
        hm.put("file", chosenFileLabel.getText());
        hm.put("paramFile", chosenParamFileLabel.getText());
        hm.put("temperature", Double.parseDouble(tempTextField.getText()));
        hm.put("field", Double.parseDouble(fieldTextField.getText()));
        hm.put("nucleus", nucTextField.getText());
        hm.put("tau", Double.parseDouble(tauTextField.getText()));
        hm.put("pressure", Double.parseDouble(pTextField.getText()));
        hm.put("mode", modeTextField.getText());
        hm.put("fitmode", fitModeChoice.getSelectionModel().getSelectedItem().toLowerCase());
        hm.put("B1field", Double.parseDouble(B1TextField.getText()));
        hm.put("Tex", Double.parseDouble(TexTextField.getText()));
        hm.put("errMode", errModeChoice.getSelectionModel().getSelectedItem());
        hm.put("errValue", Double.parseDouble(errPercentTextField.getText()));
        String[] xvals = xValTextArea.getText().split("\t");
        ArrayList<Double> fxvals = new ArrayList();
        try {
            for (int i = 0; i < xvals.length; i++) {
                fxvals.add(Double.parseDouble(xvals[i]));
            }
        } catch (NumberFormatException nfe) {
            fxvals = null;
        }
        hm.put("vcpmg", fxvals);
        dataList.add(hm);
        String fileTail = chosenFileLabel.getText().substring(0, chosenFileLabel.getText().lastIndexOf('.'));
        yamlTextField.setText(fileTail + ".yaml");
//        for (int i=0; i<dataList.size(); i++) {
//            System.out.println("dataList " + i + " " + dataList.get(i));
//        }
    }

    public void makeYAML(ActionEvent event) {
        if (dataList.isEmpty()) {
            addInfo();
        }
        makeYAML(dataList);
    }

    public void makeYAML(List data) {
        HashMap hm1 = new HashMap();
        HashMap hm2 = new HashMap();
        ArrayList<HashMap<String, Object>> dataHmList = new ArrayList();

        for (int i = 0; i < data.size(); i++) {
            HashMap hmdf = (HashMap) data.get(i);
            HashMap hmd = new HashMap(hmdf);

            String paramFile = (String) hmdf.get("paramFile");
            String paramFileName = paramFile.substring(paramFile.lastIndexOf("/") + 1, paramFile.length());
            hm2.put("mode", hmdf.get("fitmode"));
            if (!paramFileName.equals("")) {
                hm2.put("file", paramFileName);
            } else {
                hm2.put("file", "analysis.txt");
            }
//            hmdf.put("vcpmg", hmdf.get("vcpmg").toString());
            HashMap hmde = new HashMap();
            hmde.put("mode", hmdf.get("errMode"));
            hmde.put("value", hmdf.get("errValue"));
            hmdf.put("error", hmde);
//            if (hmdf.get("fitmode").equals("exp")) {
//                HashMap hmd1 = new HashMap();
//                hmd1.put("c0", 2);
//                hmd1.put("delta0", 0);
//                hmd1.put("delta", 0.008);
//                hmdf.put("delays", hmd1);
////                hmdf.remove("vcpmg");
//            }
            Set keySet = hmdf.keySet();
            if (!hmdf.get("fitmode").equals("cest")) {
                keySet.remove("Tex");
                keySet.remove("B1field");
            }
//            if (hmdf.get("fitmode").equals("cest")) {
//                keySet.remove("Tex");
//                keySet.remove("B1field");
//            }
            keySet.remove("fitmode");
            keySet.remove("paramFile");
            keySet.remove("errMode");
            keySet.remove("errValue");
            List keys = Arrays.asList(keySet);
            hmd.keySet().retainAll(keys);
            String dataFile = (String) hmdf.get("file");
            String dataFileName = dataFile.substring(dataFile.lastIndexOf("/") + 1, dataFile.length());
            hmdf.put("file", dataFileName);
            dataHmList.add(hmdf);
        }
        hm2.put("data", dataHmList);
        hm1.put("fit", hm2);

        Yaml yaml = new Yaml();
        String s = yaml.dumpAsMap(hm1);
        try (FileWriter writer = new FileWriter(yamlTextField.getText())) {
            writer.write(s);
            System.out.println(yamlTextField.getText() + " written");
        } catch (IOException ex) {
            Logger.getLogger(DataIO.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public void clearDataList(ActionEvent event) {
        dataList.clear();
    }

    public void loadInfo(ActionEvent event) {
        if (dataList.isEmpty()) {
            addInfo();
        }

        File file = null; //new File(chosenFileLabel.getText()).getAbsoluteFile();
        ResidueProperties resProp = null;
        try {
//            System.out.println(dataList);
            for (HashMap<String, Object> dataMap3 : dataList) {
                String dataFileName = (String) dataMap3.get("file");
                String fitFile = (String) dataMap3.get("paramFile");
                Double temperature = (Double) dataMap3.get("temperature");
                Double field = (Double) dataMap3.get("field");
                String nucleus = (String) dataMap3.get("nucleus");
                List<Number> vcpmgList = (List<Number>) dataMap3.get("vcpmg");
                Double tauCPMG = (Double) dataMap3.get("tau");
                tauCPMG = tauCPMG == null ? 1.0 : tauCPMG;  // fixme throw error if  ratemode and no tauCPMG
                Double B1field = (Double) dataMap3.get("B1field");
                Double Tex = (Double) dataMap3.get("Tex");
                double[] B1fieldList = new double[vcpmgList.size()];
                double[] TexList = new double[vcpmgList.size()];

                file = new File(dataFileName).getAbsoluteFile();
                Path path = file.toPath();
                Path dirPath = path.getParent();
                dataFileName = file.getName();
                String fileTail = dataFileName.substring(0, dataFileName.indexOf('.'));

//                ExperimentData expData = new ExperimentData(fileTail, nucleus, field, temperature);
                resProp = new ResidueProperties(fileTail, dataFileName);
                String expMode = (String) dataMap3.get("fitmode");
                expMode = expMode.toLowerCase();
                if (!fitFile.equals("")) {
                    resProp = DataIO.loadParametersFromFile(expMode, fitFile);
                }
                resProp.setExpMode(expMode);

                if (expMode.equals("cest")) {
                    for (int i = 0; i < B1fieldList.length; i++) {
                        B1fieldList[i] = B1field;
                        TexList[i] = Tex;
                    }
                } else {
                    B1fieldList = null;
                    TexList = null;
                }

//                System.out.println(dirPath);
//                System.out.println(dataFileName);
//                System.out.println(path);
                String textFileName = FileSystems.getDefault().getPath(dirPath.toString(), dataFileName).toString();
                String fileMode = (String) dataMap3.get("mode");
                HashMap<String, Object> errorPars = new HashMap(); //(HashMap<String, Object>) dataMap3.get("error");
                errorPars.put("mode", dataMap3.get("errMode"));
                errorPars.put("value", dataMap3.get("errValue"));
                Object delayField = dataMap3.get("delays");
                System.out.println("delays " + delayField);
                double[] delayCalc = {0.0, 0.0, 1.0};
                if (delayField instanceof Map) {
                    Map<String, Number> delayMap = (Map<String, Number>) delayField;
                    delayCalc[0] = delayMap.get("delta0").doubleValue();
                    delayCalc[1] = delayMap.get("c0").doubleValue();
                    delayCalc[2] = delayMap.get("delta").doubleValue();
                }
                System.out.println("err " + errorPars);
                if ((fileMode != null) && fileMode.equals("mpk2")) {
                    if (vcpmgList == null) {
                        ExperimentData expData = new ExperimentData(fileTail, nucleus, field, temperature, tauCPMG, null, expMode, errorPars, delayCalc, B1fieldList, TexList);
//                        loadPeakFile(textFileName, resProp, nucleus, temperature, field, tauCPMG, null, expMode, errorPars, delayCalc);
                        loadPeakFile(textFileName, expData, resProp);
                    } else {
                        double[] vcpmgs = new double[vcpmgList.size()];
                        for (int i = 0; i < vcpmgs.length; i++) {
                            vcpmgs[i] = vcpmgList.get(i).doubleValue();
                        }
                        ExperimentData expData = new ExperimentData(fileTail, nucleus, field, temperature, tauCPMG, vcpmgs, expMode, errorPars, delayCalc, B1fieldList, TexList);
//                        loadPeakFile(textFileName, resProp, nucleus, temperature, field, tauCPMG, vcpmgs, expMode, errorPars, delayCalc);
                        loadPeakFile(textFileName, expData, resProp);
                    }
                } else if (vcpmgList == null) {
                    loadTextFile(textFileName, resProp, nucleus, temperature, field);
                } else {
                    double[] vcpmgs = new double[vcpmgList.size()];
                    for (int i = 0; i < vcpmgs.length; i++) {
                        vcpmgs[i] = vcpmgList.get(i).doubleValue();
                    }
                    ExperimentData expData = new ExperimentData(fileTail, nucleus, field, temperature, tauCPMG, vcpmgs, expMode, errorPars, delayCalc, B1fieldList, TexList);
//                    loadPeakFile(textFileName, resProp, nucleus, temperature, field, tauCPMG, vcpmgs, expMode, errorPars, delayCalc);
                    loadPeakFile(textFileName, expData, resProp);
                }
            }
        } catch (IOException ioe) {
            ioe.printStackTrace();
        }
        if (file != null) {
            if (activeChart != null) {
                clearChart();
                currentResidues = null;
                simulate = false;
                fitResult = null;
            }
            XYBarChart reschartNode = PyController.mainController.getActiveChart();
            if (reschartNode == null) {
                reschartNode = PyController.mainController.addChart();

            }
//            ResidueProperties resProp = DataIO.loadParameters(fileName);
            residueProperties.put(resProp.getName(), resProp);
            String parName = "Kex";
            if (resProp.getExpMode().equals("exp")) {
                parName = "R";
            }
            ObservableList<XYChart.Series<Double, Double>> data = ChartUtil.getParMapData(resProp.getName(), "best", "0:0:0", parName);
            PyController.mainController.currentResProps = resProp;
            PyController.mainController.makeAxisMenu();
            PyController.mainController.setYAxisType(resProp.getName(), "best", "0:0:0", parName);
            reschartNode.setResProps(resProp);
            PyController.mainController.setControls();
        }
        
        fitModeChoice.setValue("Select");
        updateInfoInterface();
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
                        double aic = curveSet.getParMap().get("AIC");
                        double rms = curveSet.getParMap().get("RMS");
                        aicLabel.setText("AIC: " + aic);
                        rmsLabel.setText(" RMS: " + rms);
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
            aicLabel.setText("AIC: " + fitResult.getAicc());
            rmsLabel.setText(" RMS: " + fitResult.getRms());
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
            if (!currentStates.isEmpty()) {
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

}

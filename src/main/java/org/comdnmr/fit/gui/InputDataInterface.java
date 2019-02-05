package org.comdnmr.fit.gui;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import javafx.collections.ObservableList;
import javafx.event.ActionEvent;
import javafx.event.EventHandler;
import javafx.scene.Scene;
import javafx.scene.chart.XYChart;
import javafx.scene.control.Button;
import javafx.scene.control.CheckBox;
import javafx.scene.control.ChoiceBox;
import javafx.scene.control.Label;
import javafx.scene.control.TextArea;
import javafx.scene.control.TextField;
import javafx.scene.layout.GridPane;
import javafx.stage.DirectoryChooser;
import javafx.stage.FileChooser;
import javafx.stage.Stage;
import org.comdnmr.fit.calc.DataIO;
import org.comdnmr.fit.calc.ResidueProperties;
import static org.comdnmr.fit.gui.ChartUtil.residueProperties;
import org.controlsfx.control.textfield.TextFields;
import org.controlsfx.dialog.ExceptionDialog;
import org.yaml.snakeyaml.Yaml;

/**
 *
 * @author Martha Beckwith
 */
public class InputDataInterface {

    PyController pyController;

    GridPane inputInfoDisplay = new GridPane();
    Scene inputScene = new Scene(inputInfoDisplay, 600, 600);
    Stage infoStage = new Stage();
    TextField chosenDirLabel = new TextField();
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
    TextField yamlTextField = new TextField();
    CheckBox ppmBox = new CheckBox("ppm to Hz");
    ChoiceBox<String> errModeChoice = new ChoiceBox<>();
    TextField errPercentTextField = new TextField();
    ArrayList<HashMap<String, Object>> dataList = new ArrayList();
    Button dirChoiceButton = new Button();
    Button fileChoiceButton = new Button();
    Button xpk2ChoiceButton = new Button();
    Button paramFileChoiceButton = new Button();
    Button addButton = new Button();
    Button clearButton = new Button();
    Button yamlButton = new Button();
    Button loadButton = new Button();
    Path dirPath = null;

    public InputDataInterface(PyController controller) {
        pyController = controller;
    }

    public void inputParameters() {

        infoStage.setTitle("Input Data Parameters");
        Label fileLabel = new Label("  File:  ");
        Label dirLabel = new Label("  Directory:  ");
        Label xpk2FileLabel = new Label("  XPK2 File:  ");
        Label fitFileLabel = new Label("  Fit Parameter File:  ");
        Label fieldLabel = new Label("  B0 Field (1H MHz) :  ");
        Label tempLabel = new Label("  Temperature:  ");
        Label nucLabel = new Label("  Nucleus:  ");
        Label pLabel = new Label("  Pressure:  ");
        Label modeLabel = new Label("  Mode:  ");
        Label tauLabel = new Label("  Tau:  ");
        Label xValLabel = new Label("  X Values:  ");
        Label fitModeLabel = new Label("  Experiment Type:  ");
        Label B1FieldLabel = new Label("  B1 Field:  ");
        Label yamlLabel = new Label("  YAML File:  ");
        Label errModeLabel = new Label("  Error Mode:  ");
        Label errPercentLabel = new Label("  Error Value:  ");

        Label[] labels = {fitModeLabel, dirLabel, fileLabel, xpk2FileLabel, fitFileLabel, fieldLabel, tempLabel, nucLabel, pLabel, modeLabel,
            tauLabel, B1FieldLabel, errModeLabel, errPercentLabel, xValLabel, yamlLabel};

        dirChoiceButton.setText("Browse");
        dirChoiceButton.setOnAction(e -> chooseDirectory(e));
        chosenDirLabel.setText("");

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
        errPercentTextField.setText("5");
        errPercentTextField.setMaxWidth(textFieldWidth);
        xValTextArea.setText("");
        xValTextArea.setMaxWidth(xValAreaWidth);
        xValTextArea.setWrapText(true);
        yamlTextField.setText("");

        TextField[] texts = {fieldTextField, tempTextField, nucTextField, pTextField, modeTextField, tauTextField, B1TextField};

        inputInfoDisplay.getChildren().clear();

        for (int i = 0; i < labels.length; i++) {
            inputInfoDisplay.add(labels[i], 0, i);
        }
        for (int i = 0; i < texts.length; i++) {
            inputInfoDisplay.add(texts[i], 1, i + 5);
            texts[i].setMaxWidth(textFieldWidth);
        }

        fitModeChoice.getItems().clear();
        fitModeChoice.getItems().addAll(Arrays.asList("Select", "CPMG", "EXP", "CEST"));
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
        inputInfoDisplay.add(dirChoiceButton, 2, 1);
        inputInfoDisplay.add(chosenDirLabel, 1, 1);
        inputInfoDisplay.add(fileChoiceButton, 2, 2);
        inputInfoDisplay.add(chosenFileLabel, 1, 2);
        inputInfoDisplay.add(xpk2ChoiceButton, 2, 3);
        inputInfoDisplay.add(chosenXPK2FileLabel, 1, 3);
        inputInfoDisplay.add(paramFileChoiceButton, 2, 4);
        inputInfoDisplay.add(chosenParamFileLabel, 1, 4);
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
        yamlButton.disableProperty().bind(yamlTextField.textProperty().isEmpty());
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
        if (fitModeChoice.getSelectionModel().getSelectedItem() != null) {
            if (fitModeChoice.getSelectionModel().getSelectedItem().equals("Select")) {
                TextField[] textFields = {B1TextField, tauTextField, fieldTextField, tempTextField, nucTextField, pTextField, modeTextField,
                    errPercentTextField, yamlTextField, chosenFileLabel, chosenXPK2FileLabel, chosenParamFileLabel};
                Button[] buttons = {fileChoiceButton, xpk2ChoiceButton, paramFileChoiceButton, addButton, clearButton, loadButton};
                for (TextField textField : textFields) {
                    textField.setDisable(true);
                }
                for (Button button : buttons) {
                    button.setDisable(true);
                }
                ppmBox.setDisable(true);
                xValTextArea.setDisable(true);
                errModeChoice.setDisable(true);
            } else if (!fitModeChoice.getSelectionModel().getSelectedItem().equals("Select")) {
                TextField[] textFields = {fieldTextField, tempTextField, nucTextField, pTextField, modeTextField,
                    errPercentTextField, yamlTextField, chosenFileLabel, chosenXPK2FileLabel, chosenParamFileLabel};
                Button[] buttons = {fileChoiceButton, xpk2ChoiceButton, paramFileChoiceButton, addButton, clearButton, loadButton};
                for (TextField textField : textFields) {
                    textField.setDisable(false);
                }
                for (Button button : buttons) {
                    button.setDisable(false);
                }
                xValTextArea.setDisable(false);
                errModeChoice.setDisable(false);
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
                    tauTextField.setDisable(false);
                    ppmBox.setDisable(false);
                }
            }
        }

    }

    public void chooseDirectory(ActionEvent event) {
        DirectoryChooser dirChooser = new DirectoryChooser();
        File file = dirChooser.showDialog(infoStage);
        if (file != null) {
            chosenDirLabel.setText(file.toString());
            dirPath = file.toPath();
        }
    }

    public void chooseFile(ActionEvent event) {
        FileChooser fileChooser = new FileChooser();
        fileChooser.getExtensionFilters().addAll(new FileChooser.ExtensionFilter("mpk2 or txt File", "*.mpk2", "*.txt"));
        File file = fileChooser.showOpenDialog(infoStage);
        if (file != null) {
            Path path = dirPath.relativize(file.toPath());
            chosenFileLabel.setText(path.toString());
        }
    }

    public void chooseXPK2File(ActionEvent event) {
        FileChooser fileChooser = new FileChooser();
        fileChooser.getExtensionFilters().addAll(new FileChooser.ExtensionFilter("xpk2 File", "*.xpk2"));
        File file = fileChooser.showOpenDialog(infoStage);
        if (file != null) {
            Path path = dirPath.relativize(file.toPath());
            chosenXPK2FileLabel.setText(path.toString());

            Path path1 = file.toPath();

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
            String field = Arrays.asList(head.get(3)).get(sfInd);
            String nuc = Arrays.asList(head.get(4)).get(codeInd);
            nuc = nuc.replaceAll("[^a-zA-Z]", "");
            nucTextField.setText(nuc);
            fieldTextField.setText(field);
        }
    }

    public void chooseParamFile(ActionEvent event) {
        FileChooser fileChooser = new FileChooser();
        fileChooser.getExtensionFilters().addAll(new FileChooser.ExtensionFilter("txt File", "*.txt"));
        File file = fileChooser.showOpenDialog(infoStage);
        if (file != null) {
            String directory = file.getParent();
            String fileName = file.getName();
            chosenParamFileLabel.setText(directory + "/" + fileName);
        }
    }

    public void clearParamFile(ActionEvent event) {
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
        hm.put("B0", Double.parseDouble(fieldTextField.getText()));
        hm.put("nucleus", nucTextField.getText());
        hm.put("tau", Double.parseDouble(tauTextField.getText()));
        hm.put("pressure", Double.parseDouble(pTextField.getText()));
        hm.put("mode", modeTextField.getText());
        hm.put("fitmode", fitModeChoice.getSelectionModel().getSelectedItem().toLowerCase());
        hm.put("B1", Double.parseDouble(B1TextField.getText()));
        HashMap hmde = new HashMap();
        hmde.put("mode", errModeChoice.getSelectionModel().getSelectedItem());
        hmde.put("value", Double.parseDouble(errPercentTextField.getText()));

        hm.put("error", hmde);

        String[] xvals = xValTextArea.getText().trim().split("\t");
        if (xvals.length > 0) {
            ArrayList<Double> fxvals = new ArrayList();
            try {
                for (String xval : xvals) {
                    fxvals.add(Double.parseDouble(xval));
                }
            } catch (NumberFormatException nfe) {
                fxvals = null;
            }
            hm.put("vcpmg", fxvals);
        }
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
            Set keySet = hmd.keySet();
            if (!hmd.get("fitmode").equals("cest")) {
                keySet.remove("B1");
            }
            if ((hmd.get("vcpmg") == null) || (hmd.get("vcpmg").toString().equals(""))) {
                keySet.remove("vcpmg");
            }
            keySet.remove("fitmode");
            keySet.remove("paramFile");
            hmd.keySet().retainAll(keySet);
            String dataFile = (String) hmdf.get("file");
            String dataFileName = dataFile.substring(dataFile.lastIndexOf("/") + 1, dataFile.length());
            hmd.put("file", dataFileName);
            dataHmList.add(hmd);
        }
        hm2.put("data", dataHmList);
        hm1.put("fit", hm2);

        Yaml yaml = new Yaml();
        String s = yaml.dumpAsMap(hm1);
        Path path = FileSystems.getDefault().getPath(dirPath.toString(), yamlTextField.getText());
        try (FileWriter writer = new FileWriter(path.toFile())) {
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
        File projectDirFile = new File(chosenDirLabel.getText().trim());
        dirPath = projectDirFile.toPath();

        ResidueProperties resProp = null;
        String expMode = fitModeChoice.getSelectionModel().getSelectedItem().toLowerCase();
        resProp = new ResidueProperties(projectDirFile.getName(), projectDirFile.toString());
        expMode = expMode.toLowerCase();
        resProp.setExpMode(expMode);

        try {
            DataIO.loadDataMaps(resProp, dirPath, expMode, dataList);
        } catch (IOException ex) {
            ExceptionDialog dialog = new ExceptionDialog(ex);
            dialog.showAndWait();
            return;
        }

        if (PyController.mainController.activeChart != null) {
            PyController.mainController.clearChart();
            PyController.mainController.currentResidues = null;
            PyController.mainController.simulate = false;
            PyController.mainController.fitResult = null;
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

        fitModeChoice.setValue("Select");
        updateInfoInterface();
    }

}

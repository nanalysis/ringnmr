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
import java.util.Map;
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
import javafx.stage.FileChooser;
import javafx.stage.Stage;
import org.comdnmr.fit.calc.DataIO;
import static org.comdnmr.fit.calc.DataIO.loadPeakFile;
import static org.comdnmr.fit.calc.DataIO.loadTextFile;
import org.comdnmr.fit.calc.ExperimentData;
import org.comdnmr.fit.calc.ResidueProperties;
import static org.comdnmr.fit.gui.ChartUtil.residueProperties;
import org.controlsfx.control.textfield.TextFields;
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
    
    public InputDataInterface(PyController controller) {
        pyController = controller;
    }
    
    public void inputParameters() {

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
            TextField[] textFields = {B1TextField, tauTextField, fieldTextField, tempTextField, nucTextField, pTextField, modeTextField, TexTextField, 
                errPercentTextField, yamlTextField, chosenFileLabel, chosenXPK2FileLabel, chosenParamFileLabel};
            Button[] buttons = {fileChoiceButton, xpk2ChoiceButton, paramFileChoiceButton, addButton, clearButton, yamlButton, loadButton};
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
            TextField[] textFields = {fieldTextField, tempTextField, nucTextField, pTextField, modeTextField, TexTextField, 
                errPercentTextField, yamlTextField, chosenFileLabel, chosenXPK2FileLabel, chosenParamFileLabel};
            Button[] buttons = {fileChoiceButton, xpk2ChoiceButton, paramFileChoiceButton, addButton, clearButton, yamlButton, loadButton};
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
            for (String xval : xvals) {
                fxvals.add(Double.parseDouble(xval));
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
        }
        
        fitModeChoice.setValue("Select");
        updateInfoInterface();
    }

}

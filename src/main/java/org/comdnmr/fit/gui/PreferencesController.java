/*
 * NMRFx Processor : A Program for Processing NMR Data 
 * Copyright (C) 2004-2017 One Moon Scientific, Inc., Westfield, N.J., USA
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
package org.comdnmr.fit.gui;

import org.nmrfx.utils.properties.DirectoryOperationItem;
import org.nmrfx.utils.properties.ChoiceOperationItem;
import org.nmrfx.utils.properties.IntRangeOperationItem;
import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.ResourceBundle;
import java.util.prefs.Preferences;
import javafx.beans.value.ChangeListener;
import javafx.beans.value.ObservableValue;
import javafx.event.ActionEvent;
import javafx.fxml.FXML;
import javafx.fxml.FXMLLoader;
import javafx.fxml.Initializable;
import javafx.scene.Scene;
import javafx.scene.layout.Pane;
import javafx.stage.Stage;
import javafx.stage.StageStyle;
import org.comdnmr.fit.calc.CoMDPreferences;
import org.comdnmr.fit.calc.ExperimentData;
import org.controlsfx.control.PropertySheet;
import org.nmrfx.utils.properties.DoubleRangeOperationItem;
import org.nmrfx.utils.properties.NvFxPropertyEditorFactory;
import org.nmrfx.utils.properties.BooleanOperationItem;

/**
 *
 * @author johnsonb
 */
public class PreferencesController implements Initializable {

    @FXML
    PropertySheet prefSheet;
    ChangeListener<String> stringListener;
    ChangeListener<String> datasetListener;
    ChangeListener<String> locationListener;
    ChangeListener<Integer> nprocessListener;
    ChangeListener<Boolean> cestEqnListener;
    Stage stage;

    static File nestaNMR = null;
    static File datasetDir = null;
    private static Map<String, String> recentMap = new HashMap<>();
    static String location = null;
    static Integer nProcesses = null;
    DoubleRangeOperationItem maxFreqItem;

    @Override
    public void initialize(URL url, ResourceBundle rb) {
        prefSheet.setPropertyEditorFactory(new NvFxPropertyEditorFactory());
        prefSheet.setMode(PropertySheet.Mode.CATEGORY);
        prefSheet.setModeSwitcherVisible(false);
        prefSheet.setSearchBoxVisible(false);

        datasetListener = (ObservableValue<? extends String> observableValue, String string, String string2) -> {
            setDatasetDirectory(new File(string2.trim()));
        };
        locationListener = (ObservableValue<? extends String> observableValue, String string, String string2) -> {
            setLocation(string2.trim());
        };
        nprocessListener = (ObservableValue<? extends Integer> observableValue, Integer n1, Integer n2) -> {
            setNProcesses(n2);
        };
        cestEqnListener = (ObservableValue<? extends Boolean> observableValue, Boolean cest1, Boolean cest2) -> {
            BooleanOperationItem item = (BooleanOperationItem) observableValue;
            CoMDPreferences.setCESTEquationState(item.getName(), cest2);
        };
        ArrayList<String> locationChoices = new ArrayList<>();
        locationChoices.add("FID directory");
        locationChoices.add("Dataset directory");
        maxFreqItem = new DoubleRangeOperationItem((obs, oldV, newV) -> {
            CoMDPreferences.setCPMGMaxFreq((Double) newV);
        }, 2000.0, 100.0, 5000.0, "CPMG", "Max Freq", "Max Frequency");
        DoubleRangeOperationItem rexRatioItem = new DoubleRangeOperationItem((obs, oldV, newV) -> {
            CoMDPreferences.setRexRatio((Double) newV);
        }, 3.0, 0.0, 10.0, "CPMG", "Rex Ratio", "Rex must be this many times rmsd");
        IntRangeOperationItem nSamplesItem = new IntRangeOperationItem((obs, oldV, newV) -> {
            CoMDPreferences.setSampleSize((Integer) newV);
        }, CoMDPreferences.getSampleSize(), 10, 500, "Error Estimates", "N Samples", "Number of bootstrap samples");
        ChoiceOperationItem locationTypeItem = new ChoiceOperationItem(locationListener, getLocation(), locationChoices, "File Locations", "location", "Directory Location for Dataset");

        DirectoryOperationItem locationFileItem = new DirectoryOperationItem(datasetListener, getDatasetDirectory().getPath(), "File Locations", "Datasets", "desc");

        int nProcessesDefault = Runtime.getRuntime().availableProcessors() / 2;
        IntRangeOperationItem nProcessesItem = new IntRangeOperationItem(nprocessListener, nProcessesDefault, 1, 32, "Processor", "NProcesses", "How many parallel processes to run during processing");

        ArrayList<String> cestEqnChoices = new ArrayList<>();
        cestEqnChoices.addAll(Arrays.asList("CESTR1RHOPERTURBATIONNOEX", "CESTR1RHOPERTURBATION", "CESTR1RHOBALDWINKAY", "CESTR1RHOSD", "CESTR1RHON", "CESTR1RHOEXACT1",
                "CESTEXACT0", "CESTEXACT1", "CESTEXACT2"));
//        prefSheet.getItems().addAll(locationTypeItem, locationFileItem, nProcessesItem, maxFreqItem, rexRatioItem, nSamplesItem);
        prefSheet.getItems().addAll(nProcessesItem, maxFreqItem, rexRatioItem, nSamplesItem);
        for (String eqn : cestEqnChoices) {
            boolean defaultState = CoMDPreferences.getCESTEquationState(eqn);
            BooleanOperationItem cestEqnListItem = new BooleanOperationItem(cestEqnListener, defaultState, "CEST Equations", eqn, "List of equations to use during CEST Fitting");
            prefSheet.getItems().add(cestEqnListItem);
        }

    }

    public static PreferencesController create(Stage parent) {
        FXMLLoader loader = new FXMLLoader(PyController.class.getResource("/fxml/PreferencesScene.fxml"));
        PreferencesController controller = null;
        Stage stage = new Stage(StageStyle.DECORATED);

        try {
            Scene scene = new Scene((Pane) loader.load());
            stage.setScene(scene);
            scene.getStylesheets().add("/styles/Styles.css");

            controller = loader.<PreferencesController>getController();
            stage.setTitle("Preferences");

            stage.initOwner(parent);
            controller.stage = stage;
            stage.show();
        } catch (IOException ioE) {
            ioE.printStackTrace();
            System.out.println(ioE.getMessage());
        }

        return controller;
    }

    public Stage getStage() {
        return stage;
    }

    @FXML
    private void closeAction(ActionEvent event) {
        stage.close();
    }

    /**
     * Returns the Directory for datasets,
     *
     * @return
     */
    public static File getDatasetDirectory() {
        if (datasetDir == null) {
            Preferences prefs = Preferences.userNodeForPackage(MainApp.class);
            String filePath = prefs.get("DATASET-DIR", null);
            if (filePath != null) {
                datasetDir = new File(filePath);
            } else {
                datasetDir = new File("");
            }
        }
        return datasetDir;
    }

    /**
     * Returns the Directory for datasets,
     *
     * @param file the file or null to remove the path
     */
    public static void setDatasetDirectory(File file) {
        Preferences prefs = Preferences.userNodeForPackage(MainApp.class);
        if (file != null) {
            datasetDir = new File(file.getPath());
            prefs.put("DATASET-DIR", datasetDir.getPath());
        } else {
            datasetDir = null;
            prefs.remove("DATASET-DIR");
        }

    }

    /**
     * Returns the Directory for datasets,
     *
     * @return
     */
    public static String getLocation() {
        if (location == null) {
            Preferences prefs = Preferences.userNodeForPackage(MainApp.class);
            String value = prefs.get("LOCATION-TYPE", null);
            if (value != null) {
                location = value;
            } else {
                location = "FID Directory";
            }
        }
        return location;
    }

    /**
     * Returns the Directory for datasets,
     *
     * @param file the file or null to remove the path
     */
    public static void setLocation(String value) {
        Preferences prefs = Preferences.userNodeForPackage(MainApp.class);
        if (value != null) {
            location = value;
            prefs.put("LOCATION-TYPE", value);
        } else {
            location = null;
            prefs.remove("LOCATION-TYPE");
        }

    }

    /**
     * Saves recently opened datasets. The path is persisted in the OS specific
     * registry.
     *
     */
    public static List<Path> getRecentDatasets() {
        return getRecentFileItem("RECENT-DATASETS");
    }

    public static List<Path> getRecentProjects() {
        return getRecentFileItem("RECENT-PROJECTS");
    }

    public static void saveRecentProjects(String fileName) {
        saveRecentFileItems(fileName, "RECENT-PROJECTS");
    }

    public static void saveRecentDatasets(String fileName) {
        saveRecentFileItems(fileName, "RECENT-DATASETS");
    }

    public static void saveRecentFileItems(String fileName, String mode) {
        Preferences prefs = Preferences.userNodeForPackage(MainApp.class);
        String recentFileString = recentMap.get(mode);
        if (recentFileString == null) {
            recentFileString = prefs.get(mode, "");
            recentMap.put(mode, recentFileString);
        }
        String[] recentDatasets = recentFileString.split("\n");
        Map<String, Long> datasetMap = new HashMap<>();
        for (String recentDatasetEntry : recentDatasets) {
            String[] entry = recentDatasetEntry.split(";");
            File file = new File(entry[0]);
            if (file.exists()) {
                datasetMap.put(entry[0], Long.valueOf(entry[1]));
            }
        }
        datasetMap.put(fileName, System.currentTimeMillis());
        StringBuilder sBuilder = new StringBuilder();
        datasetMap.entrySet().stream().sorted((e1, e2) -> Long.compare(e2.getValue(), e1.getValue())).limit(15).forEach(e1 -> {
            sBuilder.append(e1.getKey());
            sBuilder.append(';');
            sBuilder.append(String.valueOf(e1.getValue()));
            sBuilder.append("\n");
        });
        recentFileString = sBuilder.toString();
        recentMap.put(mode, recentFileString);
        prefs.put(mode, recentFileString);
    }

    public static List<Path> getRecentFileItem(String mode) {
        String recentFileString = recentMap.get(mode);
        if (recentFileString == null) {
            Preferences prefs = Preferences.userNodeForPackage(MainApp.class);
            recentFileString = prefs.get(mode, "");
            recentMap.put(mode, recentFileString);
        }
        String[] recentDatasets = recentFileString.split("\n");
        List<Path> result = new ArrayList<>();
        for (String recentDatasetEntry : recentDatasets) {
            String[] entry = recentDatasetEntry.split(";");
            File file = new File(entry[0]);
            if (file.exists()) {
                Path path = file.toPath();
                result.add(path);
            }
        }
        return result;
    }

    /**
     * Returns the Directory for datasets,
     *
     * @return
     */
    public static Integer getNProcesses() {
        if (nProcesses == null) {
            Preferences prefs = Preferences.userNodeForPackage(MainApp.class);
            String value = prefs.get("NPROCESSES", null);
            if (value != null) {
                nProcesses = Integer.parseInt(value);
            } else {
                nProcesses = Runtime.getRuntime().availableProcessors() / 2;
            }
        }
        return nProcesses;
    }

    /**
     * Returns the Directory for datasets,
     *
     * @param file the file or null to remove the path
     */
    public static void setNProcesses(Integer value) {
        Preferences prefs = Preferences.userNodeForPackage(ExperimentData.class);
        if (value != null) {
            nProcesses = value;
            prefs.put("NPROCESSES", String.valueOf(value));
        } else {
            nProcesses = null;
            prefs.remove("NPROCESSES");
        }

    }

}

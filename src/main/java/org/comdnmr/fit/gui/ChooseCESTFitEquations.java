/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.fit.gui;

import java.io.IOException;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;
import java.util.ResourceBundle;
import javafx.event.ActionEvent;
import javafx.fxml.FXML;
import javafx.fxml.FXMLLoader;
import javafx.fxml.Initializable;
import javafx.geometry.Rectangle2D;
import javafx.scene.Scene;
import javafx.scene.control.Alert;
import javafx.scene.control.Button;
import javafx.scene.control.CheckBox;

import javafx.scene.layout.Pane;
import javafx.stage.Screen;
import javafx.stage.Stage;
import javafx.stage.StageStyle;
import org.comdnmr.fit.calc.CESTFit;
import org.comdnmr.fit.calc.CPMGFitResult;
import org.comdnmr.fit.calc.ResidueFitter;
import org.comdnmr.fit.calc.ResidueProperties;

/**
 * This class allows the user to choose which CEST equations to use for Fit All or Fit Selected Residues.
 *
 * @author Martha Beckwith
 *
 */
public class ChooseCESTFitEquations implements Initializable {
    
    @FXML
    private CheckBox NoExCheckBox;
    @FXML
    private CheckBox PerturbationCheckBox;
    @FXML
    private CheckBox BaldwinKayCheckBox;
    @FXML
    private CheckBox SDCheckBox;
    @FXML
    private CheckBox NCheckBox;
    @FXML
    private CheckBox R1rhoExact1CheckBox;
    @FXML
    private CheckBox Exact0CheckBox;
    @FXML
    private CheckBox Exact1CheckBox;
    @FXML
    private CheckBox Exact2CheckBox;
    @FXML
    private Button startFitsButton;
    
    Stage stage;
    CPMGFitResult fitResult;
    static boolean allRes;
    
    @Override
    public void initialize(URL url, ResourceBundle rb) {
        initializeChooser();
        MainApp.setCESTEquationController(this);
    }
    
    public static ChooseCESTFitEquations create() {
        FXMLLoader loader = new FXMLLoader(PyController.class.getResource("/fxml/CESTEquationScene.fxml"));
        ChooseCESTFitEquations controller = null;
        Stage stage = new Stage(StageStyle.DECORATED);

        try {
            Scene scene = new Scene((Pane) loader.load());
            stage.setScene(scene);
//            scene.getStylesheets().add("/styles/consolescene.css");

            controller = loader.<ChooseCESTFitEquations>getController();
            controller.stage = stage;
            stage.setTitle("CoMD/NMR CEST Equations");
            stage.show();
            Screen screen = Screen.getPrimary();
            Rectangle2D screenSize = screen.getBounds();
            stage.toFront();
            stage.setY(screenSize.getHeight() - stage.getHeight());
            ChooseCESTFitEquations consoleController = controller;
            stage.setOnCloseRequest(e -> consoleController.close());
        } catch (IOException ioE) {
            ioE.printStackTrace();
            System.out.println(ioE.getMessage());
        }

        return controller;

    }
    
    public void initializeChooser() {
        NoExCheckBox.setSelected(false);
        PerturbationCheckBox.setSelected(false);
        BaldwinKayCheckBox.setSelected(false);
        SDCheckBox.setSelected(false);
        NCheckBox.setSelected(false);
        R1rhoExact1CheckBox.setSelected(false);
        Exact0CheckBox.setSelected(false);
        Exact1CheckBox.setSelected(false);
        Exact2CheckBox.setSelected(false);
        
    }

    @FXML
    void fitEquations(ActionEvent event) {
        fitEquations(allRes);
    }
    
    void fitEquations(boolean allRes) {
        fitResult = null;
        ResidueProperties currentResProps = PyController.mainController.currentResProps;
        ResidueFitter residueFitter = PyController.mainController.residueFitter;
        List<List<String>> allResidues = new ArrayList<>();
        List<String> groupResidues = new ArrayList<>();
        ResidueChart chart = PyController.mainController.getActiveChart();
        groupResidues.addAll(chart.selectedResidues);
    
        List<String> fitEquations = new ArrayList<>();
        CheckBox[] boxes = {NoExCheckBox, PerturbationCheckBox, BaldwinKayCheckBox, SDCheckBox, NCheckBox, R1rhoExact1CheckBox, 
            Exact0CheckBox, Exact1CheckBox, Exact2CheckBox};
        String[] eqnNames = {"CESTR1RHOPERTURBATIONNOEX", "CESTR1RHOPERTURBATION", "CESTR1RHOBALDWINKAY", "CESTR1RHOSD", "CESTR1RHON", "CESTR1RHOEXACT1",
            "CESTEXACT0", "CESTEXACT1", "CESTEXACT2"};
        
        for (int i=0; i<boxes.length; i++) {
            if (boxes[i].isSelected()) {
                fitEquations.add(eqnNames[i]);
            }
        }
        
        if (fitEquations.size() > 0) {
            CESTFit.setEquationNames(fitEquations);
            PyController.mainController.equationChoice.getItems().clear();
            PyController.mainController.equationChoice.getItems().add("+");
            PyController.mainController.equationChoice.getItems().addAll(CESTFit.getEquationNames());
            PyController.mainController.equationChoice.setValue(CESTFit.getEquationNames().get(0));
            if (!allRes) {  //fit selected residues
                allResidues.add(groupResidues);
                currentResProps.setAbsValueMode(PyController.mainController.absValueModeCheckBox.isSelected());
                if (PyController.mainController.nonParBootStrapCheckBox.isSelected()) {
                    currentResProps.setBootStrapMode("nonparametric");
                } else {
                    currentResProps.setBootStrapMode("parametric");
                }
                residueFitter.fitResidues(currentResProps, allResidues);
            } else {  //fit all residues
                residueFitter.fitResidues(currentResProps);
            }
            close();
        } else {
            Alert alert = new Alert(Alert.AlertType.ERROR);
            alert.setContentText("Error: At least one fit equation must be selected.");
            alert.showAndWait();
            return;
        }
    }    
    
    public void close() {
        stage.hide();
    }

    public void show() {
        stage.show();
        stage.toFront();
    }

}

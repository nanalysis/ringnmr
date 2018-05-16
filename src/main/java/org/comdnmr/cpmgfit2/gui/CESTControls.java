/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.cpmgfit2.gui;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import javafx.fxml.FXML;
import javafx.scene.control.ChoiceBox;
import javafx.scene.control.Label;
import javafx.scene.control.Slider;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Priority;
import javafx.scene.layout.VBox;
import org.comdnmr.cpmgfit2.calc.CESTEquation;
import org.comdnmr.cpmgfit2.calc.CESTFit;
import org.comdnmr.cpmgfit2.calc.ExperimentData;
import org.comdnmr.cpmgfit2.calc.ParValueInterface;
import org.comdnmr.cpmgfit2.calc.PlotEquation;
import org.comdnmr.cpmgfit2.calc.ResidueInfo;
import org.comdnmr.cpmgfit2.calc.ResidueProperties;
import static org.comdnmr.cpmgfit2.gui.CESTControls.PARS.KEX;
import static org.comdnmr.cpmgfit2.gui.CESTControls.PARS.PB;
import static org.comdnmr.cpmgfit2.gui.CESTControls.PARS.DELTAA0;
import static org.comdnmr.cpmgfit2.gui.CESTControls.PARS.DELTAB0;
import static org.comdnmr.cpmgfit2.gui.CESTControls.PARS.R1A;
import static org.comdnmr.cpmgfit2.gui.CESTControls.PARS.R1B;
import static org.comdnmr.cpmgfit2.gui.CESTControls.PARS.R2A;
import static org.comdnmr.cpmgfit2.gui.CESTControls.PARS.R2B;
import static org.comdnmr.cpmgfit2.gui.ChartUtil.residueProperties;
import static org.comdnmr.cpmgfit2.gui.PyController.defaultField;

/**
 *
 * @author Martha Beckwith
 */
public class CESTControls implements EquationControls {

    @FXML
    ChoiceBox<String> equationSelector;
    ChoiceBox<String> stateSelector;

    String[] parNames = {"Kex", "Pb", "deltaA0", "deltaB0", "R1A", "R1B", "R2A", "R2B"};

    enum PARS implements ParControls {
        KEX("Kex", 0.0, 350.0, 50.0, 150.0),
        PB("Pb", 0.0, 1.0, 0.1, 0.1),
        DELTAA0("deltaA0", 0.0, 3000.0, 100.0, 2700.0),
        DELTAB0("deltaB0", -3000.0, 0.0, 100.0, -1250.0),
        R1A("R1A", 0.0, 4.0, 1.0, 2.5),
        R1B("R1B", 0.0, 10.0, 1.0, 5.0),
        R2A("R2A", 0.0, 40.0, 10.0, 15.0),
        R2B("R2B", 0.0, 200.0, 50.0, 120.0);

        String name;
        Slider slider;
        Label label;
        Label valueText;

        PARS(String name, double min, double max, double major, double value) {
            this.name = name;
            slider = new Slider(min, max, value);
            slider.setShowTickLabels(true);
            slider.setShowTickMarks(true);
            slider.setMajorTickUnit(major);
            label = new Label(name);
            label.setPrefWidth(60.0);
            valueText = new Label();
            valueText.setPrefWidth(60);
        }

        @Override
        public String getName() {
            return name;
        }

        @Override
        public void addTo(HBox hBox) {
            hBox.getChildren().addAll(label, slider, valueText);
            HBox.setHgrow(slider, Priority.ALWAYS);
        }

        @Override
        public Slider getSlider() {
            return slider;
        }

        @Override
        public void disabled(boolean state) {
            slider.setDisable(state);
        }

        @Override
        public void setValue(double value) {
            slider.setValue(value);
            valueText.setText(String.format("%.1f", value));
        }

        @Override
        public void setText() {
            double value = slider.getValue();
            valueText.setText(String.format("%.1f", value));
        }

        @Override
        public double getValue() {
            return slider.getValue();
        }

    }

    boolean updatingTable = false;
    PyController controller;

    @Override
    public VBox makeControls(PyController controller) {
        this.controller = controller;
        equationSelector = new ChoiceBox<>();
        equationSelector.getItems().addAll(CESTFit.getEquationNames());
        equationSelector.setValue(CESTFit.getEquationNames().get(0));
        stateSelector = new ChoiceBox<>();
        stateSelector.getItems().addAll("0:0:0", "1:0:0");
        stateSelector.setValue("0:0:0");
        VBox vBox = new VBox();
        HBox hBox1 = new HBox();
        HBox.setHgrow(hBox1, Priority.ALWAYS);
        vBox.setFillWidth(true);
        hBox1.getChildren().addAll(equationSelector, stateSelector);
        vBox.getChildren().add(hBox1);

        int i = 0;

        for (ParControls control : PARS.values()) {
            HBox hBox = new HBox();
            HBox.setHgrow(hBox, Priority.ALWAYS);
            control.addTo(hBox);

            control.getSlider().valueProperty().addListener(e -> {
                simSliderAction(control.getName());
            });
            vBox.getChildren().add(hBox);
        }

        equationSelector.valueProperty().addListener(e -> {
            equationAction();
        });
        stateSelector.valueProperty().addListener(e -> {
            stateAction();
        });
        return vBox;
    }

    void equationAction() {
        updatingTable = true;
        String equationName = equationSelector.getValue().toString();
        ResidueInfo resInfo = controller.currentResInfo;
        if (resInfo != null) {
            updatingTable = true;
            String state = stateSelector.getValue();
            List<ParValueInterface> parValues = resInfo.getParValues(equationName, state);
            controller.updateTableWithPars(parValues);
            updateSliders(parValues, equationName);
            updatingTable = false;
        }
        switch (equationName) {
            case "CESTR1RHOPERTURBATION":
                KEX.disabled(false);
                PB.disabled(false);
                DELTAA0.disabled(false);
                DELTAB0.disabled(false);
                R1A.disabled(false);
                R1B.disabled(true);
                R2A.disabled(false);
                R2B.disabled(false);
                break;
            case "CESTR1RHON":
                KEX.disabled(false);
                PB.disabled(false);
                DELTAA0.disabled(false);
                DELTAB0.disabled(false);
                R1A.disabled(false);
                R1B.disabled(true);
                R2A.disabled(false);
                R2B.disabled(true);
                break;
            case "CESTR1RHOBALDWINKAY":
                KEX.disabled(false);
                PB.disabled(false);
                DELTAA0.disabled(false);
                DELTAB0.disabled(false);
                R1A.disabled(false);
                R1B.disabled(true);
                R2A.disabled(false);
                R2B.disabled(false);
                break;
            case "CESTR1RHOSD":
                KEX.disabled(false);
                PB.disabled(false);
                DELTAA0.disabled(false);
                DELTAB0.disabled(false);
                R1A.disabled(false);
                R1B.disabled(true);
                R2A.disabled(false);
                R2B.disabled(false);
                break;
            case "CESTR1RHOEXACT1":
                KEX.disabled(false);
                PB.disabled(false);
                DELTAA0.disabled(false);
                DELTAB0.disabled(false);
                R1A.disabled(false);
                R1B.disabled(true);
                R2A.disabled(false);
                R2B.disabled(false);
                break;
            case "CESTEXACT0":
                KEX.disabled(false);
                PB.disabled(false);
                DELTAA0.disabled(false);
                DELTAB0.disabled(false);
                R1A.disabled(false);
                R1B.disabled(false);
                R2A.disabled(false);
                R2B.disabled(false);
                break;
            case "CESTEXACT1":
                KEX.disabled(false);
                PB.disabled(false);
                DELTAA0.disabled(false);
                DELTAB0.disabled(false);
                R1A.disabled(false);
                R1B.disabled(true);
                R2A.disabled(false);
                R2B.disabled(false);
                break;
            case "CESTEXACT2":
                KEX.disabled(false);
                PB.disabled(false);
                DELTAA0.disabled(false);
                DELTAB0.disabled(false);
                R1A.disabled(false);
                R1B.disabled(false);
                R2A.disabled(false);
                R2B.disabled(true);
                break;
            default:
                return;
        }
        simSliderAction(equationName);
        updatingTable = false;
    }

    void stateAction() {
        ResidueInfo resInfo = controller.currentResInfo;
        if (resInfo != null) {
            String state = stateSelector.getValue();
            if (state != null) {
                String equationName = equationSelector.getValue();
                List<ParValueInterface> parValues = resInfo.getParValues(equationName, state);
                controller.updateTableWithPars(parValues);
                updateSliders(parValues, equationName);
            }
        }

    }

    public void simSliderAction(String label) {
        if (updatingTable) {
            return;
        }
        String equationName = equationSelector.getValue().toString();
//        if (equationName.equals("CPMGSLOW") && label.equals("Rex")) {
//            return;
//        }
        double kex = KEX.getValue();
        double pb = PB.getValue();
        double deltaA0 = DELTAA0.getValue();
        double deltaB0 = DELTAB0.getValue();
        double R1a = R1A.getValue();
        double R1b = R1B.getValue();
        double R2a = R2A.getValue();
        double R2b = R2B.getValue();
        KEX.setText();
        PB.setText();
        DELTAA0.setText();
        DELTAB0.setText();
        R1A.setText();
        R1B.setText();
        R2A.setText();
        R2B.setText();
        double[] pars;
        switch (equationName) {
            case "CESTR1RHOPERTURBATION":
                pars = new double[8];
                pars[0] = kex;
                pars[1] = pb;
                pars[2] = deltaA0;
                pars[3] = deltaB0;
                pars[4] = R1a;
                pars[5] = R1b;
                pars[6] = R2a;
                pars[7] = R2b;
                break;
            case "CESTR1RHON":
                pars = new double[8];
                pars[0] = kex;
                pars[1] = pb;
                pars[2] = deltaA0;
                pars[3] = deltaB0;
                pars[4] = R1a;
                pars[5] = R1b;
                pars[6] = R2a;
                pars[7] = R2b;
                break;
            case "CESTR1RHOSD":
                pars = new double[8];
                pars[0] = kex;
                pars[1] = pb;
                pars[2] = deltaA0;
                pars[3] = deltaB0;
                pars[4] = R1a;
                pars[5] = R1b;
                pars[6] = R2a;
                pars[7] = R2b;
                break;
            case "CESTR1RHOBALDWINKAY":
                pars = new double[8];
                pars[0] = kex;
                pars[1] = pb;
                pars[2] = deltaA0;
                pars[3] = deltaB0;
                pars[4] = R1a;
                pars[5] = R1b;
                pars[6] = R2a;
                pars[7] = R2b;
                break;
            case "CESTR1RHOEXACT1":
                pars = new double[8];
                pars[0] = kex;
                pars[1] = pb;
                pars[2] = deltaA0;
                pars[3] = deltaB0;
                pars[4] = R1a;
                pars[5] = R1b;
                pars[6] = R2a;
                pars[7] = R2b;
                break;
            case "CESTEXACT0":
                pars = new double[8];
                pars[0] = kex;
                pars[1] = pb;
                pars[2] = deltaA0;
                pars[3] = deltaB0;
                pars[4] = R1a;
                pars[5] = R1b;
                pars[6] = R2a;
                pars[7] = R2b;
                break;
            case "CESTEXACT1":
                pars = new double[8];
                pars[0] = kex;
                pars[1] = pb;
                pars[2] = deltaA0;
                pars[3] = deltaB0;
                pars[4] = R1a;
                pars[5] = R1b;
                pars[6] = R2a;
                pars[7] = R2b;
                break;
            case "CESTEXACT2":
                pars = new double[8];
                pars[0] = kex;
                pars[1] = pb;
                pars[2] = deltaA0;
                pars[3] = deltaB0;
                pars[4] = R1a;
                pars[5] = R1b;
                pars[6] = R2a;
                pars[7] = R2b;
                break;
            default:
                return;
        }
//        if (equationName.equals("CPMGSLOW")) {
//            int[] map = {0, 1, 2, 3};
//            rEx = CPMGEquation.CPMGSLOW.getRex(pars, map);
//            REX.setValue(rEx);
//        }

        double[] errs = new double[pars.length];
//        int nFields = field2 > (defaultField + 10) ? 2 : 1; // addTo 10.0 to make sure slider set near to bottom gives 1 field
//        double[] fields = new double[nFields];
//        fields[0] = 1.0;
//        if (nFields > 1) {
//            fields[1] = field2 / defaultField;
//        }
        updateEquations();

        //controller.updateChartEquations(equationName, pars, errs, fields);
    }

    double[] getPars(String equationName) {
        double kex = KEX.getValue();
        double pb = PB.getValue();
        double deltaA0 = DELTAA0.getValue();
        double deltaB0 = DELTAB0.getValue();
        double R1a = R1A.getValue();
        double R1b = R1B.getValue();
        double R2a = R2A.getValue();
        double R2b = R2B.getValue();
        double[] pars;
        switch (equationName) {
            case "CESTR1RHOPERTURBATION":
                pars = new double[8];
                pars[0] = kex;
                pars[1] = pb;
                pars[2] = deltaA0;
                pars[3] = deltaB0;
                pars[4] = R1a;
                pars[5] = R1b;
                pars[6] = R2a;
                pars[7] = R2b;
                break;
            case "CESTR1RHON":
                pars = new double[8];
                pars[0] = kex;
                pars[1] = pb;
                pars[2] = deltaA0;
                pars[3] = deltaB0;
                pars[4] = R1a;
                pars[5] = R1b;
                pars[6] = R2a;
                pars[7] = R2b;
                break;
            case "CESTR1RHOSD":
                pars = new double[8];
                pars[0] = kex;
                pars[1] = pb;
                pars[2] = deltaA0;
                pars[3] = deltaB0;
                pars[4] = R1a;
                pars[5] = R1b;
                pars[6] = R2a;
                pars[7] = R2b;
                break;
            case "CESTR1RHOBALDWINKAY":
                pars = new double[8];
                pars[0] = kex;
                pars[1] = pb;
                pars[2] = deltaA0;
                pars[3] = deltaB0;
                pars[4] = R1a;
                pars[5] = R1b;
                pars[6] = R2a;
                pars[7] = R2b;
                break;
            case "CESTR1RHOEXACT1":
                pars = new double[8];
                pars[0] = kex;
                pars[1] = pb;
                pars[2] = deltaA0;
                pars[3] = deltaB0;
                pars[4] = R1a;
                pars[5] = R1b;
                pars[6] = R2a;
                pars[7] = R2b;
                break;
            case "CESTEXACT0":
                pars = new double[8];
                pars[0] = kex;
                pars[1] = pb;
                pars[2] = deltaA0;
                pars[3] = deltaB0;
                pars[4] = R1a;
                pars[5] = R1b;
                pars[6] = R2a;
                pars[7] = R2b;
                break;
            case "CESTEXACT1":
                pars = new double[8];
                pars[0] = kex;
                pars[1] = pb;
                pars[2] = deltaA0;
                pars[3] = deltaB0;
                pars[4] = R1a;
                pars[5] = R1b;
                pars[6] = R2a;
                pars[7] = R2b;
                break;
            case "CESTEXACT2":
                pars = new double[8];
                pars[0] = kex;
                pars[1] = pb;
                pars[2] = deltaA0;
                pars[3] = deltaB0;
                pars[4] = R1a;
                pars[5] = R1b;
                pars[6] = R2a;
                pars[7] = R2b;
                break;
            default:
                pars = null;
        }
        return pars;

    }

    double[] getPars(String equationName, Map<String, ParValueInterface> parValues) {
        double[] pars;
        switch (equationName) {
            case "CESTR1RHOPERTURBATION":
                pars = new double[8];
                pars[0] = parValues.get("kex").getValue();
                pars[1] = parValues.get("pb").getValue();
                pars[2] = parValues.get("deltaA0").getValue();
                pars[3] = parValues.get("deltaB0").getValue();
                pars[4] = parValues.get("R1a").getValue();
                pars[5] = parValues.get("R1b").getValue();
                pars[6] = parValues.get("R2a").getValue();
                pars[7] = parValues.get("R2b").getValue();
                break;
            case "CESTR1RHON":
                pars = new double[8];
                pars[0] = parValues.get("kex").getValue();
                pars[1] = parValues.get("pb").getValue();
                pars[2] = parValues.get("deltaA0").getValue();
                pars[3] = parValues.get("deltaB0").getValue();
                pars[4] = parValues.get("R1a").getValue();
                pars[5] = parValues.get("R1b").getValue();
                pars[6] = parValues.get("R2a").getValue();
                pars[7] = parValues.get("R2b").getValue();
                break;
            case "CESTR1RHOBALDWINKAY":
                pars = new double[8];
                pars[0] = parValues.get("kex").getValue();
                pars[1] = parValues.get("pb").getValue();
                pars[2] = parValues.get("deltaA0").getValue();
                pars[3] = parValues.get("deltaB0").getValue();
                pars[4] = parValues.get("R1a").getValue();
                pars[5] = parValues.get("R1b").getValue();
                pars[6] = parValues.get("R2a").getValue();
                pars[7] = parValues.get("R2b").getValue();
                break;
            case "CESTR1RHOSD":
                pars = new double[8];
                pars[0] = parValues.get("kex").getValue();
                pars[1] = parValues.get("pb").getValue();
                pars[2] = parValues.get("deltaA0").getValue();
                pars[3] = parValues.get("deltaB0").getValue();
                pars[4] = parValues.get("R1a").getValue();
                pars[5] = parValues.get("R1b").getValue();
                pars[6] = parValues.get("R2a").getValue();
                pars[7] = parValues.get("R2b").getValue();
                break;
            case "CESTR1RHOEXACT1":
                pars = new double[8];
                pars[0] = parValues.get("kex").getValue();
                pars[1] = parValues.get("pb").getValue();
                pars[2] = parValues.get("deltaA0").getValue();
                pars[3] = parValues.get("deltaB0").getValue();
                pars[4] = parValues.get("R1a").getValue();
                pars[5] = parValues.get("R1b").getValue();
                pars[6] = parValues.get("R2a").getValue();
                pars[7] = parValues.get("R2b").getValue();
                break;
            case "CESTEXACT0":
                pars = new double[8];
                pars[0] = parValues.get("kex").getValue();
                pars[1] = parValues.get("pb").getValue();
                pars[2] = parValues.get("deltaA0").getValue();
                pars[3] = parValues.get("deltaB0").getValue();
                pars[4] = parValues.get("R1a").getValue();
                pars[5] = parValues.get("R1b").getValue();
                pars[6] = parValues.get("R2a").getValue();
                pars[7] = parValues.get("R2b").getValue();
                break;
            case "CESTEXACT1":
                pars = new double[8];
                pars[0] = parValues.get("kex").getValue();
                pars[1] = parValues.get("pb").getValue();
                pars[2] = parValues.get("deltaA0").getValue();
                pars[3] = parValues.get("deltaB0").getValue();
                pars[4] = parValues.get("R1a").getValue();
                pars[5] = parValues.get("R1b").getValue();
                pars[6] = parValues.get("R2a").getValue();
                pars[7] = parValues.get("R2b").getValue();
                break;
            case "CESTEXACT2":
                pars = new double[8];
                pars[0] = parValues.get("kex").getValue();
                pars[1] = parValues.get("pb").getValue();
                pars[2] = parValues.get("deltaA0").getValue();
                pars[3] = parValues.get("deltaB0").getValue();
                pars[4] = parValues.get("R1a").getValue();
                pars[5] = parValues.get("R1b").getValue();
                pars[6] = parValues.get("R2a").getValue();
                pars[7] = parValues.get("R2b").getValue();
                break;
            default:
                pars = null;
        }
        return pars;

    }

    @Override
    public void updateStates(List<int[]> allStates) {
        StringBuilder sBuilder = new StringBuilder();
        stateSelector.setDisable(true);
        stateSelector.getItems().clear();
        for (int[] state : allStates) {
            sBuilder.setLength(0);
            for (int i = 1; i < state.length; i++) {
                if (sBuilder.length() > 0) {
                    sBuilder.append(':');
                }
                sBuilder.append(state[i]);
            }
            stateSelector.getItems().add(sBuilder.toString());
        }
        stateSelector.setValue(stateSelector.getItems().get(0));
        stateSelector.setDisable(false);

    }

    @Override
    public void updateSliders(List<ParValueInterface> parValues, String equationName) {
        updatingTable = true;
        for (ParValueInterface parValue : parValues) {
            String parName = parValue.getName();
            ParControls control = PARS.valueOf(parName.toUpperCase());
            if (control != null) {
                control.setValue(parValue.getValue());
            }
        }
        ResidueProperties resProps = controller.currentResProps;
//        if (resProps != null) {
//            double[] fields = resProps.getFields();
//            int iField = Integer.parseInt(stateSelector.getValue().substring(0, 1));
//            FIELD2.setValue(fields[iField]);
//            FIELD2.setText();
//        }
        equationSelector.setValue(equationName);

        updatingTable = false;
    }

    @Override
    public String getEquation() {
        return equationSelector.getValue();
    }

    void updateEquations() {
        ResidueInfo resInfo = controller.currentResInfo;
        ResidueProperties resProps = controller.currentResProps;
        List<PlotEquation> equations = new ArrayList<>();
        double[] pars;
        String equationName = equationSelector.getValue();
        //System.out.println("residueProperties = " + residueProperties);
        ResidueProperties residueProps = residueProperties.get("cest"); // fixme
        //System.out.println("expData = " + residueProps.getExperimentData("cest"));
        ExperimentData expData = residueProps.getExperimentData("cest"); // fixme
        if (resProps == null) {
            if (expData.getExtras().size() > 0) {
                pars = getPars(equationName);
                double[] errs = new double[pars.length];
                double[] extras = new double[2];
                for (int j = 0; j < expData.getExtras().size(); j++) {
                  extras[0] = 1.0;
                  extras[1] = expData.getExtras().get(j) * 2 * Math.PI;
                  //System.out.println("expData extras size = " + expData.getExtras().size()+ " extra[1] = " + extras[1]);
                  PlotEquation plotEquation = new PlotEquation(equationName, pars, errs, extras);
                  //equationCopy.setExtra(extras);

                  equations.add(plotEquation);
                }
              } else {
                  pars = getPars(equationName);
                  double[] errs = new double[pars.length];
                  double[] extras = new double[1];
                  extras[0] = 1.0;
                  PlotEquation plotEquation = new PlotEquation(equationName, pars, errs, extras);
                  //equationCopy.setExtra(extras);
                  //System.out.println("expData extras size = " + expData.getExtras().size()+ " extra[0] = " + extras[0]);
                  equations.add(plotEquation);        

              }
//            pars = getPars(equationName);
//            double[] errs = new double[pars.length];
//            double[] extras = new double[2];
//            extras[0] = 1.0;
//            extras[1] = 17.0 * 2 * Math.PI;
//            //System.out.println("updateEquations got called without resProps; extras length = "+extras.length);
//            PlotEquation plotEquation = new PlotEquation(equationName, pars, errs, extras);
//            equations.add(plotEquation);
        } else {
            double[] fields = resProps.getFields();
            String currentState = stateSelector.getValue();
            //double field2;
            if (resInfo != null) {
                for (String state : stateSelector.getItems()) {
                    int iField = Integer.parseInt(state.substring(0, 1));
                    if (state.equals(currentState)) {
                        pars = getPars(equationName);
                        //fields[iField] = FIELD2.getValue();
                    } else {
                        Map<String, ParValueInterface> parMap = new HashMap<>();
                        List<ParValueInterface> parValues = resInfo.getParValues(equationName, state);
                        for (ParValueInterface parValue : parValues) {
                            parMap.put(parValue.getName(), parValue);
                        }
                        pars = getPars(equationName, parMap);
                    }
                    if (expData.getExtras().size() > 0) {
                        pars = getPars(equationName);
                        double[] errs = new double[pars.length];
                        double[] extras = new double[2];
                        for (int j = 0; j < expData.getExtras().size(); j++) {
                          extras[0] = 1.0;
                          extras[1] = expData.getExtras().get(j) * 2 * Math.PI;
                          //System.out.println("expData extras size = " + expData.getExtras().size()+ " extra[1] = " + extras[1]);
                          PlotEquation plotEquation = new PlotEquation(equationName, pars, errs, extras);
                          //equationCopy.setExtra(extras);

                          equations.add(plotEquation);
                        }
                      } else {
                          pars = getPars(equationName);
                          double[] errs = new double[pars.length];
                          double[] extras = new double[1];
                          extras[0] = 1.0;
                          PlotEquation plotEquation = new PlotEquation(equationName, pars, errs, extras);
                          //equationCopy.setExtra(extras);
                          //System.out.println("expData extras size = " + expData.getExtras().size()+ " extra[0] = " + extras[0]);
                          equations.add(plotEquation);        

                      }
//                    double[] errs = new double[pars.length];
//                    double[] extras = new double[2];
//                    extras[0] = 1.0;
//                    extras[1] = 17.0 * 2 * Math.PI;
//                    //extras[0] = fields[iField] / fields[0];
//                    //System.out.println("updateEquations got called with resProps; extras length = "+extras.length);
//                    PlotEquation plotEquation = new PlotEquation(equationName, pars, errs, extras);
//                    equations.add(plotEquation);
                }
            }
        }
        controller.showEquations(equations);
    }
}

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
import org.comdnmr.cpmgfit2.calc.CPMGEquation;
import org.comdnmr.cpmgfit2.calc.CPMGFit;
import org.comdnmr.cpmgfit2.calc.ParValueInterface;
import org.comdnmr.cpmgfit2.calc.PlotEquation;
import org.comdnmr.cpmgfit2.calc.ResidueInfo;
import org.comdnmr.cpmgfit2.calc.ResidueProperties;
import static org.comdnmr.cpmgfit2.gui.CPMGControls.PARS.DW;
import static org.comdnmr.cpmgfit2.gui.CPMGControls.PARS.FIELD2;
import static org.comdnmr.cpmgfit2.gui.CPMGControls.PARS.KEX;
import static org.comdnmr.cpmgfit2.gui.CPMGControls.PARS.PA;
import static org.comdnmr.cpmgfit2.gui.CPMGControls.PARS.R2;
import static org.comdnmr.cpmgfit2.gui.CPMGControls.PARS.REX;
import static org.comdnmr.cpmgfit2.gui.PyController.defaultField;

/**
 *
 * @author Bruce Johnson
 */
public class CPMGControls implements EquationControls {

    @FXML
    ChoiceBox<String> equationSelector;
    ChoiceBox<String> stateSelector;

    String[] parNames = {"R2", "Kex", "Rex", "pA", "dW", "Field2"};

    enum PARS implements ParControls {
        R2("R2", 0.0, 50.0, 10.0, 10.0),
        KEX("Kex", 0.0, 4000.0, 500.0, 500.0),
        REX("Rex", 0.0, 50.0, 10.0, 8.0),
        PA("pA", 0.5, 0.99, 0.1, 0.9),
        DW("dW", 0.0, 400.0, 100.0, 100.0),
        FIELD2("Field2", 500.0, 1200.0, 100.0, 600.0);

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
            slider.setPrefWidth(200);
            label = new Label(name);
            label.setPrefWidth(50.0);
            valueText = new Label();
            valueText.setPrefWidth(50);
        }

        @Override
        public String getName() {
            return name;
        }

        @Override
        public void addTo(HBox hBox) {
            hBox.getChildren().addAll(label, slider, valueText);
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
        equationSelector.getItems().addAll(CPMGFit.getEquationNames());
        equationSelector.setValue(CPMGFit.getEquationNames().get(0));
        stateSelector = new ChoiceBox<>();
        stateSelector.getItems().addAll("0:0:0", "1:0:0");
        stateSelector.setValue("0:0:0");
        VBox vBox = new VBox();
        HBox hBox1 = new HBox();
        HBox.setHgrow(hBox1, Priority.ALWAYS);
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
            case "NOEX":
                R2.disabled(false);
                REX.disabled(true);
                KEX.disabled(true);
                PA.disabled(true);
                DW.disabled(true);
                break;
            case "CPMGFAST":
                R2.disabled(false);
                REX.disabled(false);
                KEX.disabled(false);
                PA.disabled(true);
                DW.disabled(true);
                break;
            case "CPMGSLOW":
                R2.disabled(false);
                REX.disabled(true);
                KEX.disabled(false);
                PA.disabled(false);
                DW.disabled(false);
                break;
            default:
                return;
        }
        simSliderAction(equationName);

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
        if (equationName.equals("CPMGSLOW") && label.equals("Rex")) {
            return;
        }
        double r2 = R2.getValue();
        double kEx = KEX.getValue();
        double rEx = REX.getValue();
        double pA = PA.getValue();
        double dW = DW.getValue();
        double field2 = FIELD2.getValue();
        R2.setText();
        KEX.setText();
        REX.setText();
        PA.setText();
        DW.setText();
        FIELD2.setText();
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
                break;
            default:
                return;
        }
        if (equationName.equals("CPMGSLOW")) {
            int[] map = {0, 1, 2, 3};
            rEx = CPMGEquation.CPMGSLOW.getRex(pars, map);
            REX.setValue(rEx);
        }

        double[] errs = new double[pars.length];
        int nFields = field2 > (defaultField + 10) ? 2 : 1; // addTo 10.0 to make sure slider set near to bottom gives 1 field
        double[] fields = new double[nFields];
        fields[0] = 1.0;
        if (nFields > 1) {
            fields[1] = field2 / defaultField;
        }
        updateEquations();

        //controller.updateChartEquations(equationName, pars, errs, fields);
    }

    double[] getPars(String equationName) {
        double r2 = R2.getValue();
        double kEx = KEX.getValue();
        double rEx = REX.getValue();
        double pA = PA.getValue();
        double dW = DW.getValue();
        double field2 = FIELD2.getValue();
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
                break;
            default:
                pars = null;
        }
        return pars;

    }

    double[] getPars(String equationName, Map<String, ParValueInterface> parValues) {
        double[] pars;
        switch (equationName) {
            case "NOEX":
                pars = new double[1];
                pars[0] = parValues.get("R2").getValue();
                break;
            case "CPMGFAST":
                pars = new double[3];
                pars[0] = parValues.get("Kex").getValue();
                pars[1] = parValues.get("R2").getValue();
                pars[2] = parValues.get("Rex").getValue();
                break;
            case "CPMGSLOW":
                pars = new double[4];
                pars[0] = parValues.get("Kex").getValue();
                pars[1] = parValues.get("pA").getValue();
                pars[2] = parValues.get("R2").getValue();
                pars[3] = parValues.get("dW").getValue();
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
        if (resProps != null) {
            double[] fields = resProps.getFields();
            int iField = Integer.parseInt(stateSelector.getValue().substring(0, 1));
            FIELD2.setValue(fields[iField]);
            FIELD2.setText();
        }
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
        if (resProps == null) {
            pars = getPars(equationName);
            double[] errs = new double[pars.length];
            double[] extras = new double[1];
            extras[0] = 1.0;
            PlotEquation plotEquation = new PlotEquation(equationName, pars, errs, extras);
            equations.add(plotEquation);
        } else {
            double[] fields = resProps.getFields();
            String currentState = stateSelector.getValue();
            double field2;
            if (resInfo != null) {
                for (String state : stateSelector.getItems()) {
                    int iField = Integer.parseInt(state.substring(0, 1));
                    if (state.equals(currentState)) {
                        pars = getPars(equationName);
                        fields[iField] = FIELD2.getValue();
                    } else {
                        Map<String, ParValueInterface> parMap = new HashMap<>();
                        List<ParValueInterface> parValues = resInfo.getParValues(equationName, state);
                        for (ParValueInterface parValue : parValues) {
                            parMap.put(parValue.getName(), parValue);
                        }
                        pars = getPars(equationName, parMap);
                    }
                    double[] errs = new double[pars.length];
                    double[] extras = new double[1];
                    extras[0] = fields[iField] / fields[0];
                    PlotEquation plotEquation = new PlotEquation(equationName, pars, errs, extras);
                    equations.add(plotEquation);
                }
            }
        }
        controller.showEquations(equations);
    }
}

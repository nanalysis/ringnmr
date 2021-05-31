/*
 * CoMD/NMR Software : A Program for Analyzing NMR Dynamics Data
 * Copyright (C) 2018-2019 Bruce A Johnson
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
package org.comdnmr.gui;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import javafx.fxml.FXML;
import javafx.scene.control.Label;
import javafx.scene.control.Slider;
import javafx.scene.control.TextField;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Priority;
import javafx.scene.layout.VBox;
import org.comdnmr.eqnfit.CESTEquation;
import org.comdnmr.eqnfit.CESTFitter;
import org.comdnmr.data.Experiment;
import org.comdnmr.eqnfit.ParValueInterface;
import org.comdnmr.data.ExperimentResult;
import org.comdnmr.data.ExperimentSet;
import static org.comdnmr.gui.CESTControls.PARS.KEX;
import static org.comdnmr.gui.CESTControls.PARS.PB;
import static org.comdnmr.gui.CESTControls.PARS.DELTAA0;
import static org.comdnmr.gui.CESTControls.PARS.DELTAB0;
import static org.comdnmr.gui.CESTControls.PARS.R1A;
import static org.comdnmr.gui.CESTControls.PARS.R1B;
import static org.comdnmr.gui.CESTControls.PARS.R2A;
import static org.comdnmr.gui.CESTControls.PARS.R2B;
import static org.comdnmr.gui.CESTControls.PARS.B1FIELD;
import static org.comdnmr.gui.CESTControls.PARS.TEX;
import org.comdnmr.eqnfit.CESTFitFunction;
import org.comdnmr.util.CoMDPreferences;
import javafx.scene.input.KeyCode;

/**
 *
 * @author Martha Beckwith
 */
public class CESTControls extends EquationControls {

    @FXML

    String[] parNames = {"Kex", "Pb", "deltaA0", "deltaB0", "R1A", "R1B", "R2A", "R2B", "B1field", "Tex"};

    static PyController controller = PyController.mainController;
    static final double DELTA_MIN = -20.0;
    static final double DELTA_MAX = 20.0;

    enum PARS implements ParControls {
        KEX("Kex", 0.0, 1000.0, 100.0, 150.0, 1),
        PB("Pb", 0.0, 1.0, 0.1, 0.1, 3),
        DELTAA0("deltaA0", DELTA_MIN, DELTA_MAX, 5.0, 3.0, 3),
        DELTAB0("deltaB0", DELTA_MIN, DELTA_MAX, 5.0, -2.0, 3),
        R1A("R1A", 0.0, 10.0, 1.0, 2.5, 2),
        R1B("R1B", 0.0, 10.0, 1.0, 2.5, 2),
        R2A("R2A", 0.0, 200.0, 50.0, 15.0, 2),
        R2B("R2B", 0.0, 200.0, 50.0, 120.0, 2),
        B1FIELD("B1field", 0.0, 100.0, 20.0, 20.0, 2),
        TEX("Tex", 0.0, 1.0, 0.1, 0.3, 2);

        String name;
        Slider slider;
        Label label;
        TextField valueText;
        final String format;

        PARS(String name, double min, double max, double major, double value, int digits) {
            this.name = name;
            slider = new Slider(min, max, value);
            slider.setShowTickLabels(true);
            slider.setShowTickMarks(true);
            slider.setMajorTickUnit(major);
            label = new Label(name);
            label.setPrefWidth(60.0);
            valueText = new TextField();
            valueText.setPrefWidth(60);
            format = "%." + digits + "f";

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
        public TextField getTextField() {
            return valueText;
        }

        @Override
        public void disabled(boolean state) {
            slider.setDisable(state);
        }

        @Override
        public void setValue(double value) {
            slider.setValue(value);
            setText();
        }

        @Override
        public void setText() {
            double value = slider.getValue();
            valueText.setText(String.format(format, value));
        }

        @Override
        public double getValue() {
            return slider.getValue();
        }
    }

    public List<String> getEquationNameList() {
        return CESTFitter.getEquationNames();
    }

    public void updateDeltaLimits() {
        updateDeltaLimits(DELTA_MIN, DELTA_MAX);
    }

    public void updateDeltaLimits(double min, double max) {
        PARS.DELTAA0.updateLimits(min, max);
        PARS.DELTAB0.updateLimits(min, max);
    }

    boolean updatingTable = false;
//    PyController controller;

    @Override
    public VBox makeControls(PyController controller) {
        this.controller = controller;
        VBox vBox = init();
        equationSelector.getItems().addAll(CESTEquation.getAllEquationNames());
        equationSelector.setValue(CESTEquation.getAllEquationNames()[0]);
        stateSelector.getItems().addAll("0:0:0", "1:0:0");
        stateSelector.setValue("0:0:0");
        HBox hBox1 = new HBox();
        HBox.setHgrow(hBox1, Priority.ALWAYS);
        vBox.setFillWidth(true);
        hBox1.getChildren().addAll(equationSelector, stateSelector, nucleiSelector);
        vBox.getChildren().add(hBox1);

        int i = 0;

        for (ParControls control : PARS.values()) {
            HBox hBox = new HBox();
            HBox.setHgrow(hBox, Priority.ALWAYS);
            control.addTo(hBox);

            control.getSlider().valueProperty().addListener(e -> {
                simSliderAction(control.getName());
                control.setText();
            });

            control.getTextField().setOnKeyReleased(event -> {
                if (event.getCode() == KeyCode.ENTER) {
                    String text = control.getTextField().textProperty().get();
                    if (!text.equals("")) {
                        try {
                            double value = Double.parseDouble(text);
                            control.adjustLimits(value);
                            control.getSlider().setValue(value);
                        } catch (NumberFormatException nfe) {

                        }
                    }
                }
            });

            vBox.getChildren().add(hBox);
        }

        if (controller.simulate == true) {
            equationAction();
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
        if (equationName == "") {
            equationName = equationSelector.getItems().get(0);
        }
        //System.out.println("eqnAction eqnName = " + equationName);
        ExperimentResult resInfo = controller.chartInfo.getResult();
        if (resInfo != null) {
            updatingTable = true;
            String state = stateSelector.getValue();
            List<ParValueInterface> parValues = resInfo.getParValues(equationName, state);
            controller.updateTableWithPars(parValues);
            updateSliders(parValues, equationName);
            updatingTable = false;
        }
        switch (equationName) {
            case "TROTT_PALMER":
            case "BALDWINKAY":
            case "SD":
            case "EIGENEXACT1":
            case "EXACT1":
                KEX.disabled(false);
                KEX.valueText.setDisable(false);
                PB.disabled(false);
                PB.valueText.setDisable(false);
                DELTAA0.disabled(false);
                DELTAB0.disabled(false);
                DELTAB0.valueText.setDisable(false);
                R1A.disabled(false);
                R1B.disabled(true);
                R1B.valueText.setDisable(true);
                R2A.disabled(false);
                R2B.disabled(false);
                R2B.valueText.setDisable(false);
                B1FIELD.disabled(false);
                TEX.disabled(false);
                break;
            case "LAGUERRE":
                KEX.disabled(false);
                KEX.valueText.setDisable(false);
                PB.disabled(false);
                PB.valueText.setDisable(false);
                DELTAA0.disabled(false);
                DELTAB0.disabled(false);
                DELTAB0.valueText.setDisable(false);
                R1A.disabled(false);
                R1B.disabled(true);
                R1B.valueText.setDisable(true);
                R2A.disabled(false);
                R2B.disabled(true);
                R2B.valueText.setDisable(true);
                B1FIELD.disabled(false);
                TEX.disabled(false);
                break;
            case "EXACT0":
                KEX.disabled(false);
                KEX.valueText.setDisable(false);
                PB.disabled(false);
                PB.valueText.setDisable(false);
                DELTAA0.disabled(false);
                DELTAB0.disabled(false);
                DELTAB0.valueText.setDisable(false);
                R1A.disabled(false);
                R1B.disabled(false);
                R1B.valueText.setDisable(false);
                R2A.disabled(false);
                R2B.disabled(false);
                R2B.valueText.setDisable(false);
                B1FIELD.disabled(false);
                TEX.disabled(false);
                break;
            case "EXACT2":
                KEX.disabled(false);
                KEX.valueText.setDisable(false);
                PB.disabled(false);
                PB.valueText.setDisable(false);
                DELTAA0.disabled(false);
                DELTAB0.disabled(false);
                DELTAB0.valueText.setDisable(false);
                R1A.disabled(false);
                R1B.disabled(false);
                R1B.valueText.setDisable(false);
                R2A.disabled(false);
                R2B.disabled(true);
                R2B.valueText.setDisable(true);
                B1FIELD.disabled(false);
                TEX.disabled(false);
                break;
            case "NOEX":
                KEX.disabled(true);
                KEX.valueText.setDisable(true);
                PB.disabled(true);
                PB.valueText.setDisable(true);
                DELTAA0.disabled(false);
                DELTAB0.disabled(true);
                DELTAB0.valueText.setDisable(true);
                R1A.disabled(false);
                R1B.disabled(true);
                R1B.valueText.setDisable(true);
                R2A.disabled(false);
                R2B.disabled(true);
                R2B.valueText.setDisable(true);
                B1FIELD.disabled(false);
                TEX.disabled(false);
                break;
            default:
                return;
        }
        simSliderAction(equationName);
        updatingTable = false;
    }

    void stateAction() {
        ExperimentResult resInfo = controller.chartInfo.getResult();
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
//        System.out.println("simSlider");
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
        double B1field = B1FIELD.getValue();
        double Tex = TEX.getValue();
        double[][] pars;
        switch (equationName) {
            case "TROTT_PALMER":
            case "LAGUERRE":
            case "SD":
            case "BALDWINKAY":
            case "EIGENEXACT1":
            case "EXACT0":
            case "EXACT1":
            case "EXACT2":
                pars = new double[2][8];
                pars[0][0] = kex;
                pars[0][1] = pb;
                pars[0][2] = deltaA0;
                pars[0][3] = deltaB0;
                pars[0][4] = R1a;
                pars[0][5] = R1b;
                pars[0][6] = R2a;
                pars[0][7] = R2b;
                pars[1][0] = B1field;
                pars[1][1] = Tex;
                break;
            case "NOEX":
                pars = new double[2][8];
//                pars[0][0] = kex;
//                pars[0][1] = pb;
                pars[0][0] = deltaA0;
//                pars[0][3] = deltaB0;
                pars[0][1] = R1a;
//                pars[0][5] = R1b;
                pars[0][2] = R2a;
//                pars[0][7] = R2b;
                pars[1][0] = B1field;
                pars[1][1] = Tex;
                break;
            default:
                return;
        }
//        if (equationName.equals("CPMGSLOW")) {
//            int[] map = {0, 1, 2, 3};
//            rEx = CPMGEquation.CPMGSLOW.getRex(pars, map);
//            REX.setValue(rEx);
//        }

//        double[] errs = new double[pars.length];
//        int nFields = field2 > (defaultField + 10) ? 2 : 1; // addTo 10.0 to make sure slider set near to bottom gives 1 field
//        double[] fields = new double[nFields];
//        fields[0] = 1.0;
//        if (nFields > 1) {
//            fields[1] = field2 / defaultField;
//        }
        updateEquations();

        //controller.updateChartEquations(equationName, pars, errs, fields);
    }

    public double[] getExtras() {
        double B1field = B1FIELD.getValue();
        double Tex = TEX.getValue();
        double[] extras = {B1field, Tex};
        return extras;
    }

    double[][] getPars(String equationName) {
        double kex = KEX.getValue();
        double pb = PB.getValue();
        double deltaA0 = DELTAA0.getValue();
        double deltaB0 = DELTAB0.getValue();
        double R1a = R1A.getValue();
        double R1b = R1B.getValue();
        double R2a = R2A.getValue();
        double R2b = R2B.getValue();
        double B1field = B1FIELD.getValue();
        double Tex = TEX.getValue();
        double[][] pars;
        switch (equationName) {
            case "TROTT_PALMER":
            case "LAGUERRE":
            case "SD":
            case "BALDWINKAY":
            case "EIGENEXACT1":
            case "EXACT0":
            case "EXACT1":
            case "EXACT2":
                pars = new double[2][8];
                pars[0][0] = kex;
                pars[0][1] = pb;
                pars[0][2] = deltaA0;
                pars[0][3] = deltaB0;
                pars[0][4] = R1a;
                pars[0][5] = R1b;
                pars[0][6] = R2a;
                pars[0][7] = R2b;
                pars[1][0] = B1field;
                pars[1][1] = Tex;
                break;
            case "NOEX":
                pars = new double[2][8];
//                pars[0][0] = kex;
//                pars[0][1] = pb;
                pars[0][0] = deltaA0;
//                pars[0][3] = deltaB0;
                pars[0][1] = R1a;
//                pars[0][5] = R1b;
                pars[0][2] = R2a;
//                pars[0][7] = R2b;
                pars[1][0] = B1field;
                pars[1][1] = Tex;
                break;
            default:
                pars = null;
        }
        return pars;

    }

    double[][] getPars(String equationName, Map<String, ParValueInterface> parValues) {
        double[][] pars;
        switch (equationName) {
            case "TROTT_PALMER":
            case "LAGUERRE":
            case "BALDWINKAY":
            case "SD":
            case "EIGENEXACT1":
            case "EXACT0":
            case "EXACT1":
            case "EXACT2":
                pars = new double[2][8];
                pars[0][0] = parValues.get("kex").getValue();
                pars[0][1] = parValues.get("pb").getValue();
                pars[0][2] = parValues.get("deltaA0").getValue();
                pars[0][3] = parValues.get("deltaB0").getValue();
                pars[0][4] = parValues.get("R1a").getValue();
                pars[0][5] = parValues.get("R1b").getValue();
                pars[0][6] = parValues.get("R2a").getValue();
                pars[0][7] = parValues.get("R2b").getValue();
                pars[1][0] = parValues.get("B1field").getValue();
                pars[1][1] = parValues.get("Tex").getValue();
                break;
            case "NOEX":
                pars = new double[2][8];
//                pars[0][0] = parValues.get("kex").getValue();
//                pars[0][1] = parValues.get("pb").getValue();
                pars[0][0] = parValues.get("deltaA0").getValue();
//                pars[0][3] = parValues.get("deltaB0").getValue();
                pars[0][1] = parValues.get("R1a").getValue();
//                pars[0][5] = parValues.get("R1b").getValue();
                pars[0][2] = parValues.get("R2a").getValue();
//                pars[0][7] = parValues.get("R2b").getValue();
                pars[1][0] = parValues.get("B1field").getValue();
                pars[1][1] = parValues.get("Tex").getValue();
                break;
            default:
                pars = null;
        }
        return pars;

    }

    public double[] sliderGuess(String equationName, int[][] map) {
        double kex = KEX.getValue();
        double pb = PB.getValue();
        double deltaA0 = DELTAA0.getValue();
        double deltaB0 = DELTAB0.getValue();
        double R1a = R1A.getValue();
        double R1b = R1B.getValue();
        double R2a = R2A.getValue();
        double R2b = R2B.getValue();
        double B1field = B1FIELD.getValue();
        double Tex = TEX.getValue();
        int nPars = CESTFitFunction.getNPars(map);
        double[] guesses = new double[nPars];
        switch (equationName) {
            case "TROTT_PALMER":
            case "SD":
            case "BALDWINKAY":
            case "EIGENEXACT1":
            case "EXACT1":
                for (int id = 0; id < map.length; id++) {
                    guesses[map[id][0]] = kex; //kex
                    guesses[map[id][1]] = pb; //pb
                    guesses[map[id][2]] = deltaA0; //deltaA
                    guesses[map[id][3]] = deltaB0; //deltaB
                    guesses[map[id][4]] = R1a; //R1A
                    guesses[map[id][5]] = R1a; //R1B
                    guesses[map[id][6]] = R2a; //R2A
                    guesses[map[id][7]] = R2b; //R2B
                }
                break;
            case "EXACT0":
                for (int id = 0; id < map.length; id++) {
                    guesses[map[id][0]] = kex; //kex
                    guesses[map[id][1]] = pb; //pb
                    guesses[map[id][2]] = deltaA0; //deltaA
                    guesses[map[id][3]] = deltaB0; //deltaB
                    guesses[map[id][4]] = R1a; //R1A
                    guesses[map[id][5]] = R1b; //R1B
                    guesses[map[id][6]] = R2a; //R2A
                    guesses[map[id][7]] = R2b; //R2B
                }
                break;
            case "EXACT2":
                for (int id = 0; id < map.length; id++) {
                    guesses[map[id][0]] = kex; //kex
                    guesses[map[id][1]] = pb; //pb
                    guesses[map[id][2]] = deltaA0; //deltaA
                    guesses[map[id][3]] = deltaB0; //deltaB
                    guesses[map[id][4]] = R1a; //R1A
                    guesses[map[id][5]] = R1b; //R1B
                    guesses[map[id][6]] = R2a; //R2A
                    guesses[map[id][7]] = R2a; //R2B
                }
                break;
            case "LAGUERRE":
                for (int id = 0; id < map.length; id++) {
                    guesses[map[id][0]] = kex; //kex
                    guesses[map[id][1]] = pb; //pb
                    guesses[map[id][2]] = deltaA0; //deltaA
                    guesses[map[id][3]] = deltaB0; //deltaB
                    guesses[map[id][4]] = R1a; //R1A
                    guesses[map[id][5]] = R1a; //R1B
                    guesses[map[id][6]] = R2a; //R2A
                    guesses[map[id][7]] = R2a; //R2B
                }
                break;
            case "NOEX":
                for (int id = 0; id < map.length; id++) {
//                    guesses[map[id][0]] = kex; //kex
//                    guesses[map[id][1]] = pb; //pb
                    guesses[map[id][0]] = deltaA0; //deltaA
//                    guesses[map[id][3]] = deltaA0; //deltaB
                    guesses[map[id][1]] = R1a; //R1A
//                    guesses[map[id][5]] = R1a; //R1B
                    guesses[map[id][2]] = R2a; //R2A
//                    guesses[map[id][7]] = R2a; //R2B
                }
                break;
        }
//        for(int i=0; i<guesses.length; i++){
//            System.out.println("slider guess: " + i + " " + guesses[i]);
//        }

        return guesses;

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
                if (control.getName().equals("deltaA0") || control.getName().equals("deltaB0")) {
                    if (controller.getCurrentExperimentSet() != null) {
                        if (controller.getCurrentExperimentSet().getExperimentData() != null) {
//                            double[] xVals = controller.currentResProps.getExperimentData().stream().findFirst().get().getXVals();
                            List<ExperimentResult> resVals = controller.getCurrentExperimentSet().getExperimentResults();
                            Collections.sort(resVals, (a,b) -> Integer.compare(a.getResNum(), b.getResNum()));
                            String resNum = String.valueOf(resVals.get(0).getResNum());
                            double[] xVals = controller.getCurrentExperimentSet().getExperimentData().stream().findFirst().get().getResidueData(resNum).getXValues()[0];
                            if (xVals != null) {
                                double min = Math.floor(xVals[1] / 2) * 2;
                                double max = Math.ceil(xVals[xVals.length - 1] / 2) * 2;
                                control.updateLimits(min, max);
                            }
                        }
                    }
                }
            }
        }
        ExperimentSet experimentSet = controller.getCurrentExperimentSet();
//        if (experimentSet != null) {
//            double[] fields = experimentSet.getFields();
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

    public List<String> getParNames() {
        String equationName = equationSelector.getValue().toString();
        List<String> parNames1 = new ArrayList<>();
        switch (equationName) {
            case "TROTT_PALMER":
            case "BALDWINKAY":
            case "SD":
            case "EIGENEXACT1":
            case "EXACT1":
                parNames1.add("Kex");
                parNames1.add("Pb");
                parNames1.add("deltaA0");
                parNames1.add("deltaB0");
                parNames1.add("R1A");
                parNames1.add("R2A");
                parNames1.add("R2B");
                break;
            case "EXACT0":
                parNames1.add("Kex");
                parNames1.add("Pb");
                parNames1.add("deltaA0");
                parNames1.add("deltaB0");
                parNames1.add("R1A");
                parNames1.add("R1B");
                parNames1.add("R2A");
                parNames1.add("R2B");
                break;
            case "EXACT2":
                parNames1.add("Kex");
                parNames1.add("Pb");
                parNames1.add("deltaA0");
                parNames1.add("deltaB0");
                parNames1.add("R1A");
                parNames1.add("R1B");
                parNames1.add("R2A");
                break;
            case "LAGUERRE":
                parNames1.add("Kex");
                parNames1.add("Pb");
                parNames1.add("deltaA0");
                parNames1.add("deltaB0");
                parNames1.add("R1A");
                parNames1.add("R2A");
                break;
            case "NOEX":
//                parNames1.add("Kex");
//                parNames1.add("Pb");
                parNames1.add("deltaA0");
                parNames1.add("R1A");
                parNames1.add("R2A");
                break;
        }
        return parNames1;
    }

    void updateEquations() {
//        System.out.println("CEST Controls updateEqns called.");
        ExperimentResult resInfo = controller.chartInfo.getResult();
        ExperimentSet experimentSet = controller.getCurrentExperimentSet();
        List<GUIPlotEquation> equations = new ArrayList<>();
        double[] pars;
        double[] extras1;
        String equationName = equationSelector.getValue();
        Optional<Experiment> optionalData = Optional.empty();
        if (resInfo != null) {
            if (experimentSet != null) {
                optionalData = experimentSet.getExperimentData().stream().findFirst();
                if (optionalData.isPresent() && optionalData.get().getExtras().size() > 0) {
                    for (Experiment expData : experimentSet.getExperimentData()) {
                        pars = getPars(equationName)[0];
                        List<Double> dataExtras = expData.getExtras();
                        double[] errs = new double[pars.length];
                        double[] extras = new double[3];
                        extras[0] = expData.getNucleusField();
                        extras[1] = dataExtras.get(0);
                        extras[2] = dataExtras.get(1);
//                        System.out.println("resInfo Res Num = " + resInfo.getResNum());
//                        System.out.println("extras size = " + expData.getExtras().size());
                        //System.out.println("expData extras size = " + expData.getExtras().size()+ " extra[1] = " + extras[1]);
//                        System.out.println("extras[1] = " + extras[1]);
//                        System.out.println("extras[2] = " + extras[2]);
                        GUIPlotEquation plotEquation = new GUIPlotEquation("cest", equationName, pars, errs, extras);
                        //equationCopy.setExtra(extras);

                        equations.add(plotEquation);
                    }
                } else {
                    pars = getPars(equationName)[0];
                    double[] errs = new double[pars.length];
                    double[] extras = new double[1];
                    extras[0] = CoMDPreferences.getRefField() * getNucleus().getFreqRatio(); // fixme
                    GUIPlotEquation plotEquation = new GUIPlotEquation("cest", equationName, pars, errs, extras);
                    //equationCopy.setExtra(extras);
                    //System.out.println("expData extras size = " + expData.getExtras().size()+ " extra[0] = " + extras[0]);
                    equations.add(plotEquation);
                }
            }

        } else {
            pars = getPars(equationName)[0];
            extras1 = getPars(equationName)[1];
            double[] errs = new double[pars.length];
            double[] extras = new double[3];
            extras[0] = CoMDPreferences.getRefField() * getNucleus().getFreqRatio(); // fixme
            extras[1] = extras1[0]; //17.0 * 2 * Math.PI;
            extras[2] = extras1[1]; //0.3;
            //System.out.println("updateEquations got called without experimentSet; extras length = "+extras.length);
            GUIPlotEquation plotEquation = new GUIPlotEquation("cest", equationName, pars, errs, extras);
            equations.add(plotEquation);
        }

        controller.showEquations(equations);
    }
}

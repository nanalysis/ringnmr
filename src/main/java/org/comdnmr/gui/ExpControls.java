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
import java.util.List;
import javafx.scene.control.Label;
import javafx.scene.control.Slider;
import javafx.scene.control.TextField;
import javafx.scene.input.KeyCode;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Priority;
import javafx.scene.layout.VBox;
import org.comdnmr.eqnfit.ExpFitter;
import org.comdnmr.eqnfit.ParValueInterface;
import org.comdnmr.data.ExperimentResult;
import org.comdnmr.data.ExperimentSet;
import static org.comdnmr.gui.ExpControls.PARS.A;
import static org.comdnmr.gui.ExpControls.PARS.C;
import static org.comdnmr.gui.ExpControls.PARS.R;
import org.comdnmr.eqnfit.ExpFitFunction;
import org.comdnmr.eqnfit.ExpEquation;

/**
 *
 * @author Bruce Johnson
 */
public class ExpControls extends EquationControls {

    String[] parNames = {"A", "R", "C"};

    enum PARS implements ParControls {
        A("A", 0.0, 500.0, 100.0, 100.0, 2),
        R("R", 0.0, 10.0, 2.0, 2.0, 2),
        C("C", 0.0, 100.0, 20.0, 10.0, 2),;
        String name;
        Slider slider;
        Label label;
        TextField valueText;
        String format;

        PARS(String name, double min, double max, double major, double value, int digits) {
            this.name = name;
            slider = new Slider(min, max, value);
            slider.setShowTickLabels(true);
            slider.setShowTickMarks(true);
            slider.setMajorTickUnit(major);
            label = new Label(name);
            label.setPrefWidth(50.0);
            valueText = new TextField();
            valueText.setPrefWidth(50);
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

    boolean updatingTable = false;
    PyController controller;

    public VBox makeControls(PyController controller) {
        this.controller = controller;
        VBox vBox = init();
        equationSelector.getItems().addAll(ExpEquation.getAllEquationNames());
        equationSelector.setValue(ExpEquation.getAllEquationNames()[0]);

        HBox hBox1 = new HBox();
        HBox.setHgrow(hBox1, Priority.ALWAYS);
        hBox1.getChildren().add(equationSelector);
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
        return vBox;
    }

    void equationAction() {
        String equationName = getEquation(); //equationSelector.getValue();
        switch (equationName) {
            case "EXPAB":
                A.disabled(false);
                R.disabled(false);
                C.disabled(true);
                C.valueText.setDisable(true);
                break;
            case "EXPABC":
                A.disabled(false);
                R.disabled(false);
                C.disabled(false);
                C.valueText.setDisable(false);
                break;
            default:
                return;
        }
        // simSliderAction(equationName);

    }

    public void simSliderAction(String label) {
        if (updatingTable) {
            return;
        }
        String equationName = equationSelector.getValue().toString();
        if (equationName.equals("CPMGSLOW") && label.equals("Rex")) {
            return;
        }
        updateEquations();
    }

    public double[] sliderGuess(String equationName, int[][] map) {
        double a = A.getValue();
        double r = R.getValue();
        double c = C.getValue();
        int nPars = ExpFitFunction.getNPars(map);
        double[] guesses = new double[nPars];
        switch (equationName) {
            case "EXPAB":
                for (int id = 0; id < map.length; id++) {
                    guesses[map[id][0]] = a;
                    guesses[map[id][1]] = r;
                }
                break;
            case "EXPABC":
                for (int id = 0; id < map.length; id++) {
                    guesses[map[id][0]] = a;
                    guesses[map[id][1]] = r;
                    guesses[map[id][2]] = c;
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

    }

    public void updateSliders(List<ParValueInterface> parValues, String equationName) {
        updatingTable = true;
        for (ParValueInterface parValue : parValues) {
            String parName = parValue.getName();
            ParControls control = PARS.valueOf(parName.toUpperCase());
            if (control != null) {
                control.setValue(parValue.getValue());
            }
        }
        equationSelector.setValue(equationName);
        updatingTable = false;
    }

    public String getEquation() {
        String equationName = equationSelector.getValue();
        if (equationName == "") {
            equationName = equationSelector.getItems().get(0);
        }
        return equationName; //equationSelector.getValue();
    }

    public List<String> getParNames() {
        String equationName = equationSelector.getValue().toString();
        List<String> parNames1 = new ArrayList<>();
        switch (equationName) {
            case "EXPAB":
                parNames1.add("A");
                parNames1.add("R");
                break;
            case "EXPABC":
                parNames1.add("A");
                parNames1.add("R");
                parNames1.add("C");
                break;
        }
        return parNames1;
    }

    public double[] getExtras() {
        double[] extras = {};
        return extras;
    }

    double[] getPars(String equationName) {
        double a = A.getValue();
        double r = R.getValue();
        double c = C.getValue();
        double[] pars;
        switch (equationName) {
            case "EXPAB":
                pars = new double[2];
                pars[0] = a;
                pars[1] = r;

                break;
            case "EXPABC":
                pars = new double[3];
                pars[0] = a;
                pars[1] = r;
                pars[2] = c;

                break;
            default:
                pars = null;
        }
        return pars;

    }

    void updateEquations() {
        ExperimentResult resInfo = controller.currentResInfo;
        ExperimentSet experimentSet = controller.getCurrentExperimentSet();
        List<GUIPlotEquation> equations = new ArrayList<>();
        double[] pars;
        String equationName = getEquation(); //equationSelector.getValue();
        if (experimentSet == null) {
            pars = getPars(equationName);
            double[] errs = new double[pars.length];
            double[] extras = new double[1];
            extras[0] = 1.0;
            GUIPlotEquation plotEquation = new GUIPlotEquation("exp", equationName, pars, errs, extras);
            equations.add(plotEquation);
        } else {
            pars = getPars(equationName);
            double[] errs = new double[pars.length];
            double[] extras = new double[1];
            extras[0] = 1.0;
            GUIPlotEquation plotEquation = new GUIPlotEquation("exp", equationName, pars, errs, extras);
            equations.add(plotEquation);
        }
        controller.showEquations(equations);
    }
}

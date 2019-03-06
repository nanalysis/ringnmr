/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.fit.gui;

import java.util.ArrayList;
import java.util.List;
import javafx.scene.control.Label;
import javafx.scene.control.Slider;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Priority;
import javafx.scene.layout.VBox;
import org.comdnmr.fit.calc.ExpFit;
import org.comdnmr.fit.calc.ParValueInterface;
import org.comdnmr.fit.calc.PlotEquation;
import org.comdnmr.fit.calc.ResidueInfo;
import org.comdnmr.fit.calc.ResidueProperties;
import static org.comdnmr.fit.gui.ExpControls.PARS.A;
import static org.comdnmr.fit.gui.ExpControls.PARS.C;
import static org.comdnmr.fit.gui.ExpControls.PARS.R;
import org.comdnmr.fit.calc.CalcExpDecay;

/**
 *
 * @author Bruce Johnson
 */
public class ExpControls extends EquationControls {

    String[] parNames = {"A", "R", "C"};

    enum PARS implements ParControls {
        A("A", 0.0, 500.0, 100.0, 100.0),
        R("R", 0.0, 10.0, 2.0, 2.0),
        C("C", 0.0, 100.0, 20.0, 10.0),;
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
        
        @Override
        public void updateLimits(double min, double max) {
            slider.setMin(min);
            slider.setMax(max);
        }

    }

    boolean updatingTable = false;
    PyController controller;

    public VBox makeControls(PyController controller) {
        this.controller = controller;
        VBox vBox = init();
        equationSelector.getItems().addAll(ExpFit.getEquationNames());
        equationSelector.setValue(ExpFit.getEquationNames().get(0));
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
                break;
            case "EXPABC":
                A.disabled(false);
                R.disabled(false);
                C.disabled(false);
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
        A.setText();
        R.setText();
        C.setText();
        updateEquations();
    }

    public double[] sliderGuess(String equationName, int[][] map) {
        double a = A.getValue();
        double r = R.getValue();
        double c = C.getValue();
        int nPars = CalcExpDecay.getNPars(map);
        double[] guesses = new double[nPars];
        switch (equationName) {
            case "EXPAB":
                for (int id = 0; id < map.length; id++) {
                    guesses[map[id][0]] = a;
                    guesses[map[id][1]] = r;
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
                pars[2] = r;

                break;
            default:
                pars = null;
        }
        return pars;

    }

    void updateEquations() {
        ResidueInfo resInfo = controller.currentResInfo;
        ResidueProperties resProps = controller.currentResProps;
        List<GUIPlotEquation> equations = new ArrayList<>();
        double[] pars;
        String equationName = getEquation(); //equationSelector.getValue();
        if (resProps == null) {
            pars = getPars(equationName);
            double[] errs = new double[pars.length];
            double[] extras = new double[1];
            extras[0] = 1.0;
            GUIPlotEquation plotEquation = new GUIPlotEquation(equationName, pars, errs, extras);
            equations.add(plotEquation);
        } else {
            pars = getPars(equationName);
            double[] errs = new double[pars.length];
            double[] extras = new double[1];
            extras[0] = 1.0;
            GUIPlotEquation plotEquation = new GUIPlotEquation(equationName, pars, errs, extras);
            equations.add(plotEquation);
        }
        controller.showEquations(equations);
    }
}

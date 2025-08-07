package org.comdnmr.gui;

import java.util.ArrayList;
import java.util.List;
import java.util.Optional;

import org.comdnmr.data.Experiment;
import org.comdnmr.data.ExperimentResult;
import org.comdnmr.data.ExperimentSet;
import org.comdnmr.eqnfit.ParValueInterface;
import org.comdnmr.eqnfit.SSR1RhoEquation;
import org.comdnmr.util.CoMDPreferences;

import javafx.scene.control.Label;
import javafx.scene.control.Slider;
import javafx.scene.control.TextField;
import javafx.scene.input.KeyCode;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Priority;
import javafx.scene.layout.VBox;

public class SSR1RhoControls extends EquationControls {

    String[] parNames = {"tauc", "s2", "omega1", "omegaR"};

    static PyController controller = PyController.mainController;
    static final double OMEGA_1_MIN = 2.0 * Math.PI * 2.0e3;  // ν1 (min) = 2kHz
    static final double OMEGA_1_MAX = 2.0 * Math.PI * 40.0e3;  // ν1 (max) = 40kHz
    static final double OMEGA_R_MIN = 2.0 * Math.PI * 2.0e3;  // νR (min) = 2kHz
    static final double OMEGA_R_MAX = 2.0 * Math.PI * 40.0e3;  // νR (max) = 40kHz

    // TODO this is largely parroting the relevant code in R1RhoControls...
    // Should this be generalized?
    enum PARS implements ParControls {
        TAUC("τc", 0.0, 1.0e-4, 1.0e-6, 1.0e-5, 10),
        S2("S²", 0.0, 1.0, 0.01, 0.3285, 3),
        // TODO currently, ω1 is fixed as the independent variable
        // OMEGA1("ω1", OMEGA_1_MIN, OMEGA_1_MAX, 1.0e3, 2.0 * Math.PI * 10.0e3, 2),
        OMEGAR("ωR", OMEGA_R_MIN, OMEGA_R_MAX, 1.0e3, 2.0 * Math.PI * 10.0e3, 2);

        String name;
        Slider slider;
        Label label;
        TextField valueText;
        final String format;

        // >>>>>>>>>>>>>>>
        // TODO: Identical to R1RhoControls...
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
        // <<<<<<<<<<<<<<<
    }

    boolean updatingTable = false;

    // >>>>>>>>>>>>>>>>
    // TODO: This is all the same as R1RhoControls, except for the use of
    // SSR1RhoEquation.getAllEquationNames
    public VBox makeControls(PyController controller) {
        // TODO: Warning:
        // The static field SSR1RhoControls.controller should be accessed in a
        // static way (Java 570425420)
        this.controller = controller;
        // Calls `EquationControls.init()`
        VBox vBox = init();
        String[] eqNames = SSR1RhoEquation.getAllEquationNames();
        equationSelector.getItems().addAll(eqNames);
        equationSelector.setValue(eqNames[0]);

        // TODO: verify this is correct
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

            control.getTextField().setOnKeyReleased(e -> {
                if (e.getCode() == KeyCode.ENTER) {
                    String text = control.getTextField().textProperty().get();
                    if (!text.equals("")) {
                        try {
                            double value = Double.parseDouble(text);
                            control.adjustLimits(value);
                            control.getSlider().setValue(value);
                        } catch (NumberFormatException exc) {
                            // Do Nothing
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
    // <<<<<<<<<<<<<<<<

    public double[] getExtras() {
        double omegaR = PARS.OMEGAR.getValue();
        double[] extras = {omegaR};
        return extras;
    }

    void equationAction() {
        // >>>>>>>>>>
        // TODO: All the same as R1RhoControls...
        String equationName = equationSelector.getValue().toString();
        if (equationName == "") {
            equationName = equationSelector.getItems().get(0);
        }
        ExperimentResult resInfo = controller.chartInfo.getResult();
        updatingTable = true;
        if (resInfo != null) {
            String state = stateSelector.getValue();
            List<ParValueInterface> parValues = resInfo.getParValues(equationName, state);
            controller.updateTableWithPars(parValues);
            updateSliders(parValues, equationName);
            updatingTable = false;
        }
        // <<<<<<<<<<

    switch (equationName) {
            case "CSA":
                for (PARS par: PARS.values()) {
                    par.disabled(false);
                    par.valueText.setDisable(false);
                }
                break;

            default:
                return;
        }

        simSliderAction(equationName);
        updatingTable = false;
    }

    // >>>>>>>>>>>
    // TODO: Identical to R1RhoControls...
    void stateAction() {
        ExperimentResult resInfo = controller.chartInfo.getResult();
        if (resInfo != null) {
            String state = stateSelector.getValue();
            if (state != null) {
                List<ParValueInterface> parValues = resInfo.getParValues(getEquation(), state);
                controller.updateTableWithPars(parValues);
                updateSliders(parValues, getEquation());
            }
        }
    }
    // <<<<<<<<<<<

    @Override
    String getEquation() {
        return equationSelector.getValue();
    }

    // TODO: Why is this public, while other methods are package-private?
    public void simSliderAction(String label) {
        if (updatingTable) {
            return;
        }

        updateEquations();
    }

    void updateEquations() {
        ExperimentResult eRes = controller.chartInfo.getResult();
        ExperimentSet eSet = controller.getCurrentExperimentSet();
        List<GUIPlotEquation> equations = new ArrayList<>();
        Optional<Experiment> optionalData = Optional.empty();

        double[] pars;
        double[] extras;
        double[] errs;

        if (eRes != null && eSet != null) {
            // TODO: I believe his is applicable to real data
            // See line 687 onwards in R1RhoControls
        } else {
            double[][] allPars = getPars(getEquation());
            pars = allPars[0];
            double[] extras1 = allPars[1];
            errs = new double[pars.length];
            extras = new double[2];
            extras[0] = extras1[0]; // omegaR
            extras[1] =  CoMDPreferences.getRefField() * getNucleus().getFreqRatio();

            GUIPlotEquation plotEquation = new GUIPlotEquation("ssr1rho", getEquation(), pars, errs, extras);
            equations.add(plotEquation);
        }

        controller.showEquations(equations);
    }

    double[][] getPars(String equationName) {
        double[][] pars;
        switch (equationName) {
            case "CSA":
                double tauc = PARS.TAUC.getValue();
                double s2 = PARS.S2.getValue();
                double omegaR = PARS.OMEGAR.getValue();
                pars = new double[2][2];
                pars[0][0] = tauc;
                pars[0][1] = s2;
                pars[1][0] = omegaR;
                break;
            default:
                pars = null;
                break;
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
            stateSelector.getItems().add(stateSelector.getItems().get(0));
            stateSelector.setDisable(false);
        }
    }

    public double[] sliderGuess(String equationName, int[][] map) {
        double tauc = PARS.TAUC.getValue();
        double s2 = PARS.S2.getValue();
        double[] guesses = new double[2];
        switch (equationName) {
            case "CSA":
                for (int id = 0; id < map.length; id++) {
                    guesses[map[id][0]] = tauc;
                    guesses[map[id][1]] = s2;
                }
                break;
        }
        return guesses;
    }

    public List<String> getParNames() {
        return List.of(parNames[0], parNames[1]);
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
        equationSelector.setValue(equationName);
        updatingTable = false;
    }
}

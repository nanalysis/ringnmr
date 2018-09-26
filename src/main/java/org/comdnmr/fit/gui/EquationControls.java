/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.fit.gui;

import java.util.List;
import javafx.scene.layout.VBox;
import org.comdnmr.fit.calc.ParValueInterface;

/**
 *
 * @author Bruce Johnson
 */
public interface EquationControls {

    String getEquation();

    List<String> getParNames();

    public double[] getExtras();

    VBox makeControls(PyController controller);

    void updateSliders(List<ParValueInterface> parValues, String equationName);

    void updateStates(List<int[]> allStates);

    public void simSliderAction(String label);

    public double[] sliderGuess(String equationName, int[][] map);

}

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.cpmgfit2.gui;

import java.util.List;
import javafx.scene.layout.VBox;
import org.comdnmr.cpmgfit2.calc.ParValueInterface;

/**
 *
 * @author Bruce Johnson
 */
public interface EquationControls {

    String getEquation();

    VBox makeControls(PyController controller);

    void updateSliders(List<ParValueInterface> parValues, String equationName);

    public void simSliderAction(String label);

}

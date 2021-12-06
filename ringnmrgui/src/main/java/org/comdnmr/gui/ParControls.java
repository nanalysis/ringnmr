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

import javafx.scene.control.Slider;
import javafx.scene.control.TextField;
import javafx.scene.layout.HBox;

/**
 *
 * @author Bruce Johnson
 */
public interface ParControls {

    String getName();

    void addTo(HBox hBox);

    Slider getSlider();

    TextField getTextField();

    void disabled(boolean state);

    void setValue(double value);

    void setText();

    double getValue();

    default void updateLimits(double min, double max) {
        Slider slider = getSlider();
        slider.setMin(min);
        slider.setMax(max);
        double range = max - min;
        double incValue = range / 5;
        slider.setMajorTickUnit(incValue);
    }

    default void adjustLimits(double value) {
        Slider slider = getSlider();
        double curMax = slider.getMax();
        double curMin = slider.getMin();
        if (value > curMax) {
            double newMax = curMax * 2.0;
            while (newMax < value) {
                newMax *= 2.0;
            }
            updateLimits(curMin, newMax);
            slider.setValue(value);
        } else if (value < curMax / 5.0) {
            double newMax = curMax / 5.0;
            updateLimits(curMin, newMax);
            slider.setValue(value);
        }

    }

}

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

import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;
import javafx.scene.paint.Color;
import javafx.stage.FileChooser;

/**
 *
 * @author brucejohnson
 */
public class SecondaryStructure {

    int start;
    int end;
    int mode;
    String label;
    Color color = Color.GRAY;

    public SecondaryStructure(int start, int end, int mode, String label) {
        this.start = start;
        this.end = end;
        this.mode = mode;
        this.label = label;
        if (this.label.equals(".")) {
            this.label = "";
        }
    }

    public SecondaryStructure(int start, int end, String strMode, String label, String colorName) {
        this.start = start;
        this.end = end;
        if (strMode.startsWith("h") || strMode.startsWith("H")) {
            this.mode = 1;
        } else if (strMode.startsWith("s") || strMode.startsWith("S")) {
            this.mode = 2;
        } else if (strMode.startsWith("d") || strMode.startsWith("D")) {
            this.mode = 3;
        } else {
            this.mode = 0;
        }
        this.label = label;
        if (this.label.equals(".")) {
            this.label = "";
        }
        this.color = Color.web(colorName);
    }

    boolean isHelix() {
        return mode == 1;
    }

    boolean isSheet() {
        return mode == 2;
    }

    boolean isCoil() {
        return mode == 0;
    }

    boolean isDisulfide() {
        return mode == 3;
    }

    int getStart() {
        return start;
    }

    int getEnd() {
        return end;
    }

    Color getColor() {
        return color;
    }

    String getLabel() {
        return label;
    }

    public static List<SecondaryStructure> loadFromFile() {
        FileChooser fileChooser = new FileChooser();
        fileChooser.getExtensionFilters().addAll(new FileChooser.ExtensionFilter("Secondary Structure File", "*.txt"));
        File file = fileChooser.showOpenDialog(null);
        List<SecondaryStructure> ssValues = new ArrayList<>();
        if (file != null) {
            try {
                List<String> lines = Files.readAllLines(file.toPath(), Charset.defaultCharset());
                for (String line : lines) {
                    line = line.trim();
                    if ((line.length() != 0) && (line.charAt(0) != '#')) {
                        String[] fields = line.split(" +");
                        if (fields.length > 2) {
                            try {
                                int start = Integer.parseInt(fields[0]);
                                int end = Integer.parseInt(fields[1]);
                                String type = fields[2];
                                String label = "";
                                String color = "gray";
                                if (fields.length > 3) {
                                    label = fields[3];
                                }
                                if (fields.length > 4) {
                                    color = fields[4];
                                }
                                SecondaryStructure ss = new SecondaryStructure(start, end, type, label, color);
                                ssValues.add(ss);
                            } catch (NumberFormatException nfE) {
                                nfE.printStackTrace();
                            }
                        }

                    }

                }
            } catch (IOException ex) {
            }
        }
        return ssValues;

    }
}

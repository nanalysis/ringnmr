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

import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import javafx.scene.canvas.Canvas;
import javafx.scene.canvas.GraphicsContext;
import javafx.scene.paint.Color;
import org.nmrfx.graphicsio.GraphicsContextInterface;
import org.nmrfx.graphicsio.GraphicsContextProxy;
import org.nmrfx.graphicsio.GraphicsIOException;

/**
 *
 * @author brucejohnson
 */
public class SSPainter {

    Canvas canvas;
    double cHeight = 40.0;
    final double height = 10;
    final double arrowEdge = 5;
    int nDisulfides = 0;

    List<SecondaryStructure> secondaryStructures;

    public SSPainter(Canvas canvas, List<SecondaryStructure> secondaryStructures) {
        this.canvas = canvas;
        this.secondaryStructures = secondaryStructures;
        int n = 0;
        for (SecondaryStructure ss : secondaryStructures) {
            if (ss.isDisulfide()) {
                n++;
            }
        }
        this.nDisulfides = n;
        double diSulfideOffset = 0.0;
        if (nDisulfides > 0) {
            diSulfideOffset = nDisulfides * arrowEdge + 3;
        }

        cHeight = 3.0 * height + diSulfideOffset;

    }

    public double getHeight() {
        return cHeight;
    }

    void paintSS(double xPos, double yPos, double cWidth, double start, double end) {
        GraphicsContext gCC = canvas.getGraphicsContext2D();
        GraphicsContextInterface gC = new GraphicsContextProxy(gCC);
        try {
            gC.clearRect(xPos, yPos, cWidth, cHeight);
            paintSS(gC, xPos, yPos, cWidth, start, end);
        } catch (GraphicsIOException ex) {
        }
    }

    void paintSS(GraphicsContextInterface gC, double xPos, double yPos, double cWidth, double resStart, double resEnd) {
        try {
            final double xScale = cWidth / (resEnd - resStart);
            //  final double center = cHeight - 1.0 * height - arrowEdge;
            final double center = yPos + height + arrowEdge;
            int nDisulfide = 0;
            for (SecondaryStructure ss : secondaryStructures) {
                Color color = ss.getColor();
                Color darkColor = color.darker();
                String label = ss.getLabel();
                double start = (ss.getStart() - resStart) * xScale - xScale / 2 + xPos;
                double end = (ss.getEnd() - resStart) * xScale + xScale / 2 + xPos;
                double labelCenterX = (start + end) / 2.0;
                double labelCenterY = center - height;
                System.out.println(center + " " + start + " " + end + " " + ss.mode);
                if (ss.isCoil()) {
                    gC.setFill(color);
                    gC.setStroke(color);
                    gC.strokeLine(start, center, end, center);
                } else if (ss.isDisulfide()) {
                    gC.setFill(color);
                    gC.setStroke(color);
                    double ssStart = start + xScale / 2;
                    double ssEnd = end - xScale / 2;
                    double offset = nDisulfide * arrowEdge + center + height + 3;
                    gC.strokeLine(ssStart, center, ssStart, offset);
                    gC.strokeLine(ssEnd, center, ssEnd, offset);
                    gC.strokeLine(ssStart, offset, ssEnd, offset);
                    nDisulfide++;
                } else if (ss.isSheet()) {
                    gC.setFill(color);
                    gC.setStroke(color);
                    gC.beginPath();
                    gC.moveTo(start, center - height / 2);
                    gC.lineTo(end - xScale, center - height / 2);
                    gC.lineTo(end - xScale, center - height / 2 - arrowEdge);
                    gC.lineTo(end, center);
                    gC.lineTo(end - xScale, center + height / 2 + arrowEdge);
                    gC.lineTo(end - xScale, center + height / 2);
                    gC.lineTo(start, center + height / 2);
                    gC.closePath();
                    gC.fill();
                } else if (ss.isHelix()) {
                    for (int hElem = 0; hElem < 3; hElem++) {
                        for (int iHelix = ss.getStart() + hElem; iHelix < ss.getEnd(); iHelix += 3) {
                            int helixOffset = iHelix - (int) resStart;
                            double x1 = helixOffset * xScale - xScale / 2 + xPos;
                            double x2 = helixOffset * xScale + xScale / 2 + xPos;
                            double y1 = center + height;
                            double y2 = center - height;
                            double dX = (x2 - x1);
                            double slope2 = 1.5 * dX / 2.0;

                            if ((hElem % 3) == 0) {
                                double x1p = x1;
                                double x2p = x2;
                                double y3p = (y1 + y2) / 2;

                                double x3p = x1p + slope2;
                                double x4p = x3p + dX;

                                gC.beginPath();
                                gC.moveTo(x1p, y1);
                                gC.lineTo(x3p, y3p);
                                gC.lineTo(x4p, y3p);
                                gC.lineTo(x2p, y1);
                                gC.closePath();
                                gC.setFill(color);
                                gC.setStroke(color);
                                gC.fill();
                                gC.stroke();
                            } else if ((hElem % 3) == 1) {
                                double x1p = x1;
                                double x2p = x2;
                                double x3p = x1p + slope2 - dX;
                                double x4p = x3p + slope2;
                                double x5p = x4p + dX;
                                double x6p = x5p - slope2;

                                double y3p = (y1 + y2) / 2;

                                gC.beginPath();
                                gC.moveTo(x3p, y3p);
                                gC.lineTo(x4p, y2);
                                gC.lineTo(x5p, y2);
                                gC.lineTo(x6p, y3p);

                                gC.closePath();
                                gC.setFill(color);
                                gC.setStroke(color);
                                gC.fill();
                                gC.stroke();

                                double x7p = x5p + slope2;
                                double x8p = x7p - dX;
                                gC.beginPath();
                                gC.moveTo(x4p, y2);
                                gC.lineTo(x5p, y2);
                                gC.lineTo(x7p, y3p);
                                gC.lineTo(x8p, y3p);

                                gC.closePath();
                                gC.setFill(darkColor);
                                gC.setStroke(darkColor);
                                gC.fill();
                                gC.stroke();

                            } else if ((hElem % 3) == 2) {

                                double x1p = x1;
                                double x2p = x2;

                                double x3p = x2p - slope2;
                                double x4p = x3p + dX;
                                double x5p = x4p + slope2;
                                double x6p = x5p - dX;
                                double y3p = (y1 + y2) / 2;

                                gC.beginPath();
                                gC.moveTo(x3p, y3p);
                                gC.lineTo(x4p, y3p);
                                gC.lineTo(x5p, y1);
                                gC.lineTo(x6p, y1);
                                gC.closePath();
                                gC.setFill(darkColor);
                                gC.setStroke(darkColor);
                                gC.fill();
                                gC.stroke();
                            }
                        }
                    }
                }
                if (!label.equals("")) {
                    // Text text = new Text(labelCenterX, labelCenterY, label);
                }
            }
        } catch (GraphicsIOException ioE) {

        }
    }

    public double getPrefHeight() {
        int nDisulfide = 0;
        for (SecondaryStructure ss : secondaryStructures) {
            if (ss.isDisulfide()) {
                nDisulfide++;

            }
        }
        double height = 10.0;
        double arrowEdge = 5.0;
        return 2.5 * height + (1 + nDisulfide) * arrowEdge;

    }

}

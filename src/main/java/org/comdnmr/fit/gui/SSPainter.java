/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.fit.gui;

import java.util.ArrayList;
import javafx.collections.ObservableList;
import javafx.scene.Group;
import javafx.scene.canvas.Canvas;
import javafx.scene.canvas.GraphicsContext;
import javafx.scene.chart.NumberAxis;
import javafx.scene.paint.Color;
import javafx.scene.shape.Line;
import javafx.scene.shape.Polygon;
import javafx.scene.text.Text;

/**
 *
 * @author brucejohnson
 */
public class SSPainter {

    Canvas canvas;
    NumberAxis xAxis;
    NumberAxis yAxis;
    ArrayList<SecondaryStructure> secondaryStructures;

    public SSPainter(Canvas canvas, NumberAxis xAxis, NumberAxis yAxis, ArrayList<SecondaryStructure> secondaryStructures) {
        this.canvas = canvas;
        this.xAxis = xAxis;
        this.yAxis = yAxis;
        this.secondaryStructures = secondaryStructures;

    }

    public SSPainter(NumberAxis xAxis, NumberAxis yAxis, ArrayList<SecondaryStructure> secondaryStructures) {
        this.xAxis = xAxis;
        this.yAxis = yAxis;
        this.secondaryStructures = secondaryStructures;

    }

    void paintSS() {
        GraphicsContext gC = canvas.getGraphicsContext2D();
        final double cWidth = xAxis.getWidth();
        final double cHeight = yAxis.getHeight();
        canvas.setWidth(cWidth);
        canvas.setHeight(cHeight);
        gC.clearRect(0, 0, cWidth, cHeight);
        final double height = 10;
        final double xScale = cWidth / xAxis.getUpperBound();
        final double arrowEdge = 5;
        final double center = cHeight - 1.0 * height - arrowEdge;
        for (SecondaryStructure ss : secondaryStructures) {
            Color color = ss.getColor();
            Color darkColor = color.darker();
            String label = ss.getLabel();
            double start = ss.getStart() * xScale - xScale / 2;
            double end = ss.getEnd() * xScale + xScale / 2;
            double labelCenterX = (start + end) / 2.0;
            double labelCenterY = center - height;

            if (ss.isSheet()) {
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
                        double x1 = iHelix * xScale - xScale / 2;
                        double x2 = iHelix * xScale + xScale / 2;
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
                Text text = new Text(labelCenterX, labelCenterY, label);
            }
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

    void addSS(Group group) {
        final double cWidth = xAxis.getWidth();
        final double xScale = cWidth / (xAxis.getUpperBound() - xAxis.getLowerBound());
        final double height = 10;
        final double arrowEdge = 5;
//        final double center = cHeight - 1.0 * height - arrowEdge;
        final double center = height;
        int firstRes = (int) xAxis.getLowerBound();

        group.getChildren().clear();
        double startX = xAxis.getLayoutX() + 5;  // 5 is a little fudge factor that helps with appearance
        double disx = xAxis.getDisplayPosition(xAxis.getLowerBound());
        double disy = yAxis.getLayoutX() + yAxis.getWidth();
        int nDisulfide = 0;
        for (SecondaryStructure ss : secondaryStructures) {
            Color color = ss.getColor();
            Color darkColor = color.darker();
            String label = ss.getLabel();
            double start = (ss.getStart() - firstRes) * xScale - xScale / 2 + startX;
            double end = (ss.getEnd() - firstRes) * xScale + xScale / 2 + startX;
            double labelCenterX = (end + start) / 2.0;
            double labelCenterY = center - height;

            if (ss.isCoil()) {
                Line line = new Line(start, center, end, center);
                group.getChildren().add(line);
            } else if (ss.isDisulfide()) {
                double ssStart = start + xScale / 2;
                double ssEnd = end - xScale / 2;
                double offset = nDisulfide * arrowEdge + center + height + 3;
                Line line = new Line(ssStart, center, ssStart, offset);
                line.setStroke(color);
                group.getChildren().add(line);
                line = new Line(ssEnd, center, ssEnd, offset);
                line.setStroke(color);
                group.getChildren().add(line);
                line = new Line(ssStart, offset, ssEnd, offset);
                line.setStroke(color);
                group.getChildren().add(line);
                nDisulfide++;
            } else if (ss.isSheet()) {
                Polygon polygon = new Polygon();
                ObservableList<Double> points = polygon.getPoints();
                points.add(start);
                points.add(center - height / 2);
                points.add(end - xScale);
                points.add(center - height / 2);
                points.add(end - xScale);
                points.add(center - height / 2 - arrowEdge);
                points.add(end);
                points.add(center);
                points.add(end - xScale);
                points.add(center + height / 2 + arrowEdge);
                points.add(end - xScale);
                points.add(center + height / 2);
                points.add(start);
                points.add(center + height / 2);
                polygon.setFill(color);
                polygon.setStroke(color);
                group.getChildren().add(polygon);

            } else if (ss.isHelix()) {
                for (int hElem = 0; hElem < 3; hElem++) {
                    for (int iHelix = ss.getStart() + hElem; iHelix < ss.getEnd(); iHelix += 3) {
                        double x1 = (iHelix - firstRes) * xScale - xScale / 2 + startX;
                        double x2 = (iHelix - firstRes) * xScale + xScale / 2 + startX;
                        double y1 = center + height;
                        double y2 = center - height;
                        double dX = (x2 - x1);
                        double slope2 = 1.5 * dX / 2.0;

                        switch (hElem % 3) {
                            case 0: {
                                double x1p = x1;
                                double x2p = x2;
                                double y3p = (y1 + y2) / 2;
                                double x3p = x1p + slope2;
                                double x4p = x3p + dX;
                                Polygon polygon = new Polygon();
                                ObservableList<Double> points = polygon.getPoints();
                                points.add(x1p);
                                points.add(y1);
                                points.add(x3p);
                                points.add(y3p);
                                points.add(x4p);
                                points.add(y3p);
                                points.add(x2p);
                                points.add(y1);
                                polygon.setFill(darkColor);
                                polygon.setStroke(darkColor);
                                group.getChildren().add(polygon);
                                break;
                            }
                            case 1: {
                                double x1p = x1;
                                double x2p = x2;
                                double x3p = x1p + slope2 - dX;
                                double x4p = x3p + slope2;
                                double x5p = x4p + dX;
                                double x6p = x5p - slope2;
                                double y3p = (y1 + y2) / 2;
                                Polygon polygon = new Polygon();
                                ObservableList<Double> points = polygon.getPoints();
                                points.add(x3p);
                                points.add(y3p);
                                points.add(x4p);
                                points.add(y2);
                                points.add(x5p);
                                points.add(y2);
                                points.add(x6p);
                                points.add(y3p);
                                polygon.setFill(darkColor);
                                polygon.setStroke(darkColor);
                                group.getChildren().add(polygon);
                                double x7p = x5p + slope2;
                                double x8p = x7p - dX;
                                Polygon polygon2 = new Polygon();
                                ObservableList<Double> points2 = polygon2.getPoints();
                                points2.add(x4p);
                                points2.add(y2);
                                points2.add(x5p);
                                points2.add(y2);
                                points2.add(x7p);
                                points2.add(y3p);
                                points2.add(x8p);
                                points2.add(y3p);
                                polygon2.setFill(color);
                                polygon2.setStroke(color);
                                group.getChildren().add(polygon2);
                                break;
                            }
                            case 2: {
                                double x1p = x1;
                                double x2p = x2;
                                double x3p = x2p - slope2;
                                double x4p = x3p + dX;
                                double x5p = x4p + slope2;
                                double x6p = x5p - dX;
                                double y3p = (y1 + y2) / 2;
                                Polygon polygon = new Polygon();
                                ObservableList<Double> points = polygon.getPoints();
                                points.add(x3p);
                                points.add(y3p);
                                points.add(x4p);
                                points.add(y3p);
                                points.add(x5p);
                                points.add(y1);
                                points.add(x6p);
                                points.add(y1);
                                polygon.setFill(color);
                                polygon.setStroke(color);
                                group.getChildren().add(polygon);
                                break;
                            }
                            default:
                                break;
                        }
                    }
                }
            }
            if (!label.equals("")) {
                Text text = new Text(label);
                text.setLayoutX(labelCenterX - text.prefWidth(-1) / 2);
                text.setLayoutY(labelCenterY);
                text.setFill(Color.BLACK);
                group.getChildren().add(text);
            }

        }
    }

}

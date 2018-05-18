package org.comdnmr.fit.gui;

/**
 * Data extra values for storing value, lowPercentile and highPercentile.
 */
public class ErrorExtraValues {

    private double value;
    private double lowPercentile;
    private double highPercentile;

    public ErrorExtraValues(double value, double lowPercentile, double highPercentile) {
        this.value = value;
        this.lowPercentile = lowPercentile;
        this.highPercentile = highPercentile;
    }

    public double getValue() {
        return value;
    }

    public double getLowPercentile() {
        return lowPercentile;
    }

    public double getHighPercentile() {
        return highPercentile;
    }

    @Override
    public String toString() {
        return value + " " + lowPercentile + " " + highPercentile;
    }
}

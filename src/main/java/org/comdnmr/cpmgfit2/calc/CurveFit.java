package org.comdnmr.cpmgfit2.calc;

import java.util.Map;

/**
 *
 * @author Bruce Johnson
 */
public class CurveFit {

    final String state;
    final String resNum;
    final Map<String, Double> parMap;
    final PlotEquation plotEquation;

    public CurveFit(String state, String resNum, Map<String, Double> parMap, PlotEquation plotEquation) {
        this.state = state;
        this.resNum = resNum;
        this.parMap = parMap;
        this.plotEquation = plotEquation;
    }

    public Map<String, Double> getParMap() {
        return parMap;
    }

    public PlotEquation getEquation() {
        return plotEquation;
    }

    public String getResNum() {
        return resNum;
    }
    
    public String getState() {
        return state;
    }
}

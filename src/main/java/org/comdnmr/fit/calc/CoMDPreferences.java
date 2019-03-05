/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.fit.calc;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.prefs.Preferences;

/**
 *
 * @author brucejohnson
 */
public class CoMDPreferences {

    static private Double cpmgMaxFreq = null;
    static private Double rexRatio = null;
    static private Double startRadius = null;
    static private Double finalRadius = null;
    static private Double tolerance = null;
    static private Boolean weightFit = null;
    static private Boolean absValueFit = null;
    static private Boolean nonParametricBootstrap = null;
    static private Integer sampleSize = null;
    static private String optimizer = null;
    static private String bootStrapOptimizer = null;
    private static Map<String, Boolean> cestEqnMap = null;
    private static Map<String, Boolean> r1rhoEqnMap = null;
    private static final String DEFAULT_CEST_EQNS = "CESTR1RHOPERTURBATIONNOEX;true\nCESTR1RHOPERTURBATION;true";
    private static final String DEFAULT_R1RHO_EQNS = "R1RHOPERTURBATIONNOEX;true\nR1RHOPERTURBATION;true";

    static Preferences getPrefs() {
        Preferences prefs = Preferences.userNodeForPackage(CoMDPreferences.class);
        return prefs;
    }

    public static Double getCPMGMaxFreq() {
        if (cpmgMaxFreq == null) {
            String value = getPrefs().get("CPMG_MAX_FREQ", "2000.0");
            cpmgMaxFreq = Double.parseDouble(value);
        }
        return cpmgMaxFreq;
    }

    public static void setCPMGMaxFreq(Double value) {
        cpmgMaxFreq = value;
        if (value != null) {
            getPrefs().put("CPMG_MAX_FREQ", value.toString());
        } else {
            getPrefs().remove("CPMG_MAX_FREQ");
        }
    }

    public static Integer getSampleSize() {
        if (sampleSize == null) {
            String value = getPrefs().get("SAMPLE_SIZE", "50");
            sampleSize = Integer.parseInt(value);
        }
        return sampleSize;
    }

    public static void setSampleSize(Integer value) {
        sampleSize = value;
        if (value != null) {
            getPrefs().put("SAMPLE_SIZE", value.toString());
        } else {
            getPrefs().remove("SAMPLE_SIZE");
        }
    }

    public static String getOptimizer() {
        if (optimizer == null) {
            String value = getPrefs().get("OPTIMIZER", "CMA-ES");
            optimizer = value;
        }
        return optimizer;
    }

    public static void setOptimizer(String value) {
        optimizer = value;
        if (value != null) {
            getPrefs().put("OPTIMIZER", value);
        } else {
            getPrefs().remove("OPTIMIZER");
        }
    }

    public static String getBootStrapOptimizer() {
        if (bootStrapOptimizer == null) {
            String value = getPrefs().get("BOOTSTRAP-OPTIMIZER", "BOBYQA");
            bootStrapOptimizer = value;
        }
        return bootStrapOptimizer;
    }

    public static void setBootStrapOptimizer(String value) {
        bootStrapOptimizer = value;
        if (value != null) {
            getPrefs().put("BOOTSTRAP-OPTIMIZER", value);
        } else {
            getPrefs().remove("BOOTSTRAP-OPTIMIZER");
        }
    }

    public static Double getRexRatio() {
        if (rexRatio == null) {
            String value = getPrefs().get("REX_RATIO", "3.0");
            rexRatio = Double.parseDouble(value);
        }
        return rexRatio;
    }

    public static void setRexRatio(Double value) {
        rexRatio = value;
        if (value != null) {
            getPrefs().put("REX_RATIO", value.toString());
        } else {
            getPrefs().remove("REX_RATIO");
        }
    }

    public static Double getStartingRadius() {
        if (startRadius == null) {
            String value = getPrefs().get("START_RADIUS", "20.0");
            startRadius = Double.parseDouble(value);
        }
        return startRadius;
    }

    public static void setStartingRadius(Double value) {
        startRadius = value;
        if (value != null) {
            getPrefs().put("START_RADIUS", value.toString());
        } else {
            getPrefs().remove("START_RADIUS");
        }
    }

    public static Double getFinalRadius() {
        if (finalRadius == null) {
            String value = getPrefs().get("FINAL_RADIUS", "-5");
            finalRadius = Double.parseDouble(value);
        }
        return finalRadius;
    }

    public static void setFinalRadius(Double value) {
        finalRadius = value;
        if (value != null) {
            getPrefs().put("FINAL_RADIUS", value.toString());
        } else {
            getPrefs().remove("FINAL_RADIUS");
        }
    }

    public static Double getTolerance() {
        if (tolerance == null) {
            String value = getPrefs().get("TOLERANCE", "-5");
            tolerance = Double.parseDouble(value);
        }
        return tolerance;
    }

    public static void setTolerance(Double value) {
        tolerance = value;
        if (value != null) {
            getPrefs().put("TOLERANCE", value.toString());
        } else {
            getPrefs().remove("TOLERANCE");
        }
    }

    public static Boolean getWeightFit() {
        if (weightFit == null) {
            String value = getPrefs().get("WEIGHT_FIT", "true");
            weightFit = Boolean.parseBoolean(value);
        }
        return weightFit;
    }

    public static void setWeightFit(Boolean value) {
        weightFit = value;
        if (value != null) {
            getPrefs().put("WEIGHT_FIT", value.toString());
        } else {
            getPrefs().remove("WEIGHT_FIT");
        }
    }

    public static Boolean getAbsValueFit() {
        if (absValueFit == null) {
            String value = getPrefs().get("ABSVAL_FIT", "true");
            absValueFit = Boolean.parseBoolean(value);
        }
        return absValueFit;
    }

    public static void setAbsValueFit(Boolean value) {
        absValueFit = value;
        if (value != null) {
            getPrefs().put("ABSVAL_FIT", value.toString());
        } else {
            getPrefs().remove("ABSVAL_FIT");
        }
    }

    public static Boolean getNonParametric() {
        if (nonParametricBootstrap == null) {
            String value = getPrefs().get("NONPARAMETRIC", "true");
            nonParametricBootstrap = Boolean.parseBoolean(value);
        }
        return nonParametricBootstrap;
    }

    public static void setNonParametric(Boolean value) {
        nonParametricBootstrap = value;
        if (value != null) {
            getPrefs().put("NONPARAMETRIC", value.toString());
        } else {
            getPrefs().remove("NONPARAMETRIC");
        }
    }

    static Map<String, Boolean> getCESTEqnMap() {
        if (cestEqnMap == null) {
            Preferences prefs = Preferences.userNodeForPackage(ExperimentData.class);
            String eqns = prefs.get("CEST_EQNS", DEFAULT_CEST_EQNS);
            cestEqnMap = stringToMap(eqns);
        }
        return cestEqnMap;
    }

    static Map<String, Boolean> getR1RhoEqnMap() {
        if (r1rhoEqnMap == null) {
            Preferences prefs = Preferences.userNodeForPackage(ExperimentData.class);
            String eqns = prefs.get("R1RHO_EQNS", DEFAULT_R1RHO_EQNS);
            r1rhoEqnMap = stringToMap(eqns);
        }
        return r1rhoEqnMap;
    }
    
    public static void setCESTEqnMap(String eqns) {
        cestEqnMap = stringToMap(eqns);
    }

    public static void setR1RhoEqnMap(String eqns) {
        r1rhoEqnMap = stringToMap(eqns);
    }

    static Map<String, Boolean> stringToMap(String s) {
        Map<String, Boolean> map = new HashMap<>();
        String[] stringList = s.split("\n");
        for (String strValue : stringList) {
            String[] strValueParts = strValue.split(";");
            if (strValueParts[1].equals("true")) {
                map.put(strValueParts[0], Boolean.TRUE);
            } else {
                map.put(strValueParts[0], Boolean.FALSE);
            }
        }
        return map;
    }

    static String mapToString(Map<String, Boolean> map) {
        StringBuilder sBuilder = new StringBuilder();
        for (String eqn : map.keySet()) {
            sBuilder.append(eqn);
            sBuilder.append(';');
            sBuilder.append(map.get(eqn));
            sBuilder.append("\n");
        }
        return sBuilder.toString();
    }

    public static void saveCESTEqnPrefs() {
        Preferences prefs = Preferences.userNodeForPackage(ExperimentData.class);
        String eqnString = DEFAULT_CEST_EQNS;
        if (cestEqnMap != null) {
            eqnString = mapToString(getCESTEqnMap());
        }
        prefs.put("CEST_EQNS", eqnString);
    }

    public static void saveR1RhoEqnPrefs() {
        Preferences prefs = Preferences.userNodeForPackage(ExperimentData.class);
        String eqnString = DEFAULT_R1RHO_EQNS;
        if (r1rhoEqnMap != null) {
            eqnString = mapToString(getR1RhoEqnMap());
        }
        prefs.put("R1RHO_EQNS", eqnString);
    }

    public static List<String> getActiveCESTEquations() {
        Map<String, Boolean> map = getCESTEqnMap();
        List<String> cestEqnList = new ArrayList<>();
        for (String eqn : map.keySet()) {
            if (map.get(eqn)) {
                cestEqnList.add(eqn);
            }
        }
        return cestEqnList;
    }

    public static List<String> getActiveR1RhoEquations() {
        Map<String, Boolean> map = getR1RhoEqnMap();
        List<String> r1rhoEqnList = new ArrayList<>();
        for (String eqn : map.keySet()) {
            if (map.get(eqn)) {
                r1rhoEqnList.add(eqn);
            }
        }
        return r1rhoEqnList;
    }

    public static void setCESTEquationState(String equation, boolean state) {
        Map<String, Boolean> map = getCESTEqnMap();
        map.put(equation, state);
        saveCESTEqnPrefs();
    }

    public static boolean getCESTEquationState(String equation) {
        Map<String, Boolean> map = getCESTEqnMap();
        boolean state = false;
        if (map.containsKey(equation) && map.get(equation)) {
            state = true;
        }
        return state;
    }

    public static void setR1RhoEquationState(String equation, boolean state) {
        Map<String, Boolean> map = getR1RhoEqnMap();
        map.put(equation, state);
        saveR1RhoEqnPrefs();
    }

    public static boolean getR1RhoEquationState(String equation) {
        Map<String, Boolean> map = getR1RhoEqnMap();
        boolean state = false;
        if (map.containsKey(equation) && map.get(equation)) {
            state = true;
        }
        return state;
    }

}

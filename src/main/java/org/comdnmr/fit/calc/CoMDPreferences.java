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
    static private Integer sampleSize = null;
    private static Map<String, Boolean> eqnMap = new HashMap<>();

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
    
    public static void setEqnMap(Map<String, Boolean> map) {
        eqnMap = map;
    }
    public static void saveEqnPrefs() {
        Preferences prefs = Preferences.userNodeForPackage(ExperimentData.class);
        StringBuilder sBuilder = new StringBuilder();
        for (String eqn : eqnMap.keySet()){
            sBuilder.append(eqn);
            sBuilder.append(';');
            sBuilder.append(eqnMap.get(eqn));
            sBuilder.append("\n");
        }
        prefs.put("CESTEQNS", sBuilder.toString());
    }

    public static List<String> getEqns() {
        Preferences prefs = Preferences.userNodeForPackage(ExperimentData.class);
        String eqns = prefs.get("CESTEQNS", null);
        String[] cestEqnPrefList = eqns.split("\n");
        List<String> cestEqnList = new ArrayList<>();
        for (String cestEqnPref : cestEqnPrefList) {
            String[] cestEqnPrefs = cestEqnPref.split(";");
            if (cestEqnPrefs[1].equals("true")) {
                cestEqnList.add(cestEqnPrefs[0]);
            }
        }
        return cestEqnList;
    }

}

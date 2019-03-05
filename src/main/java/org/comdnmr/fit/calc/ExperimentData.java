package org.comdnmr.fit.calc;

import java.util.HashMap;
import java.util.Set;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author Bruce Johnson
 */
public class ExperimentData {

    public enum Nuclei {

        H1("H", 1, "1/2", 99.98, 1.00000) {
        },
        F19("F", 19, "1/2", 100, 0.94077) {
        },
        P31("P", 31, "1/2", 100, 0.40481) {
        },
        C13("C", 13, "1/2", 1.108, 0.25144) {
        },
        N15("N", 15, "1/2", 0.37, 0.10133) {
        };
        String name;
        int num;
        int spin;
        double abundance;
        double freqRatio;

        public double getRatio() {
            return freqRatio;
        }

        Nuclei(final String name, final int num, final String spin, final double abundance, final double freqRatio) {
            this.name = name;
            this.num = num;
            this.spin = Integer.parseInt(spin.substring(0, 1));
            this.abundance = abundance;
            this.freqRatio = freqRatio;
        }

        public static Nuclei get(String name) {
            Nuclei nuc = null;
            try {
                nuc = valueOf(name);
            } catch (IllegalArgumentException iaE) {
                for (Nuclei testNuc : values()) {
                    if (testNuc.name.equals(name)) {
                        nuc = testNuc;
                        break;
                    }
                }
            }
            return nuc;
        }

    }

    String name;
    HashMap<String, ResidueData> residueData;
    final double field;
    final double nucleusField;
    final double temperature;
    final Double tau;
    final double[] xvals;
    final String expMode;
    final HashMap<String, Object> errorPars;
    final double[] delayCalc;
    final Double B1field;
    double errFraction = 0.05;
    final String nucleusName;
    List<Double> extras = new ArrayList<>();
    private String state = "";

    public ExperimentData(String name, String nucleus, double field, double temperature, Double tau, double[] xvals,
            String expMode, HashMap<String, Object> errorPars, double[] delayCalc, Double B1field) {
        this.name = name;
        this.field = field;
        this.temperature = temperature;
        this.tau = tau;
        this.xvals = xvals;
        this.expMode = expMode;
        this.errorPars = errorPars;
        this.delayCalc = delayCalc;
        this.B1field = B1field;
        if (nucleus == null) {
            nucleus = "H1";
        }
        this.nucleusName = nucleus;
        double ratio = 1.0;
        Nuclei nuc = Nuclei.get(nucleusName);
        if (nuc != null) {
            ratio = nuc.getRatio();
        }
        nucleusField = field * ratio;
        residueData = new HashMap<>();
    }

    @Override
    public String toString() {
        StringBuilder sBuilder = new StringBuilder();
        sBuilder.append(name).append(" ").append(field).append(" ").
                append(B1field).append(" nres ").append(residueData.size());
        return sBuilder.toString();
    }

    public String getName() {
        return name;
    }

    public void setErrFraction(double errFraction) {
        this.errFraction = errFraction;
    }

    public double getErrFraction() {
        return (errFraction);
    }

    public void addResidueData(String resNum, ResidueData data) {
        residueData.put(resNum, data);
    }

    public ResidueData getResidueData(String resNum) {
        return residueData.get(resNum);
    }

    public Set<String> getResidues() {
        return residueData.keySet();
    }

    public double getField() {
        return field;
    }

    public double getNucleusField() {
        return nucleusField;
    }

    /**
     * @return the state
     */
    public String getState() {
        return state;
    }

    /**
     * @param state the state to set
     */
    public void setState(String state) {
        this.state = state;
    }

    public void setExtras(List extravals) {
        this.extras.addAll(extravals);
    }

    public List<Double> getExtras() {
        return extras;
    }

    public Double getTau() {
        return tau;
    }

    public double[] getXVals() {
        return xvals;
    }

    public String getExpMode() {
        return expMode;
    }

    public HashMap<String, Object> getErrorPars() {
        return errorPars;
    }

    public double[] getDelayCalc() {
        return delayCalc;
    }

    public Double getB1Field() {
        return B1field;
    }
}

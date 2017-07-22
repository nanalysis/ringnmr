package org.comdnmr.cpmgfit2.calc;

import java.util.Arrays;
import java.util.ArrayList;
import java.util.DoubleSummaryStatistics;
import java.util.stream.IntStream;
import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.SimpleBounds;
import org.apache.commons.math3.optim.SimpleValueChecker;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.CMAESOptimizer;
import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.SynchronizedRandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.util.FastMath;

public class CalcRDisp implements MultivariateFunction {

    long startTime = 0;
    double[] xValues;
    double[] fieldValues;
    double[] yValues;
    double[] errValues;
    double[] fields;
    int[] idNums;
    int[][] map;
    int nID = 1;
    int reportAt = 10;
    CPMGEquation equation;
    boolean absMode = false;
    int nSim = 200;
    boolean reportFitness = false;
    double[] rexErrors = new double[nID];
    static RandomGenerator random = new SynchronizedRandomGenerator(new Well19937c());

    public static enum CPMGEquation {
        NOEX("noex", 0, "R2") {
            @Override
            public double calculate(double[] par, int[] map, double x, int idNum, double field) {
                /*
                 Zero exchange cpmg function
                 R2(1/tcp)=R2
                 */

                double R2 = par[map[0]];
                double value = R2;
                return value;
            }

            @Override
            double[] guess(double[] xValues, double[] yValues, int[] idNums, int nID, double field) {
                double[] guesses = new double[nID];
                for (int id = 0; id < nID; id++) {
                    double mean = DataUtil.getMeanValue(yValues, idNums, id);
                    guesses[id] = mean;
                }
                return guesses;
            }

            @Override
            double[][] boundaries(double[] xValues, double[] yValues, int[] idNums, int nID, double field) {
                double[] guesses = guess(xValues, yValues, idNums, nID, field);
                double[][] boundaries = new double[2][guesses.length];

                for (int id = 0; id < nID; id++) {
                    boundaries[0][id] = 0.0;
                    boundaries[1][id] = guesses[id] * 4;
                }
                return boundaries;
            }

            @Override
            public double getRex(double[] pars) {
                return 0.0;
            }

            @Override
            double getKex(double[] pars) {
                return 0.0;
            }

            @Override
            double getRex(double[] pars, int id) {
                return 0.0;
            }

            @Override
            double getKex(double[] pars, int id) {
                return 0.0;
            }

            @Override
            public int[][] makeMap(int n) {
                int[][] map = new int[n][1];
                for (int i = 0; i < n; i++) {
                    map[i][0] = i;
                }
                return map;
            }

            @Override
            public int[][] makeMap(int n, int m) {
                int[][] map = new int[n][1];
                for (int i = 0; i < n; i++) {
                    map[i][0] = i;
                }
                return map;
            }
        },
        CPMGFAST("cpmgfast", 1, "Kex", "R2", "Rex") {

            @Override
            public double calculate(double[] par, int[] map, double vu, int idNum, double field) {
                /*
                 Fast limit cpmg function
                 R2(1/tcp)=R2+Rex*(1 - (2*Tau/tcp)*tanh(1/(2*Tau/tcp)))
                 */
// do we need a different R2 for each field
                double kEx = par[map[0]];
                double R2 = par[map[1]];
                double Rex = par[map[2]];
                double value;
                if (kEx <= 0.0) {
                    value = R2;
                } else {
                    // rex = pa *pb *dw *dw / kEx;
                    double tauCP = 1.0 / (2.0 * vu);
                    double fieldAdjust = field / fieldRef;
                    Rex *= fieldAdjust * fieldAdjust;
                    value = R2 + Rex * (1 - 2.0 * FastMath.tanh(0.5 * kEx * tauCP) / (kEx * tauCP));
                }
                return value;
            }

            @Override
            double[] guess(double[] xValues, double[] yValues, int[] idNums, int nID, double field) {
                double[] guesses = new double[1 + nID * 2];
                double kExSum = 0.0;
                for (int id = 0; id < nID; id++) {
                    double minY = DataUtil.getMinValue(yValues, idNums, id);
                    double maxY = DataUtil.getMaxValue(yValues, idNums, id);

                    double r2 = minY;
                    double rex = maxY - minY;
                    if (rex < 0.0) {
                        rex = 0.0;
                    }
                    guesses[1 + id * 2] = r2;
                    guesses[2 + id * 2] = rex;
                    double vMid = DataUtil.getMidValue(yValues, xValues, idNums, id);
                    double tauMid = 1.0 / (2.0 * vMid);
                    if (rex >= 0) {
                        kExSum += 1.915 / (0.5 * tauMid); // 1.915 comes from solving equation iteratively at tcp rex 0.5 half max
                    }
                }
                guesses[0] = kExSum /= nID;
                return guesses;

            }

            @Override
            double[][] boundaries(double[] xValues, double[] yValues, int[] idNums, int nID, double field) {
                double[] guesses = guess(xValues, yValues, idNums, nID, field);
                double[][] boundaries = new double[2][guesses.length];
                boundaries[0][0] = 0.0;
                boundaries[1][0] = guesses[0] * 4;
                for (int id = 0; id < nID; id++) {
                    boundaries[0][id * 2 + 1] = 0.0;
                    boundaries[1][id * 2 + 1] = guesses[id * 2 + 1] * 4;
                    boundaries[0][id * 2 + 2] = 0.0;
                    boundaries[1][id * 2 + 2] = guesses[id * 2 + 2] * 4;

                }
                return boundaries;
            }

            @Override
            public double getRex(double[] pars) {
                return pars[2];
            }

            @Override
            double getRex(double[] pars, int id) {
                return pars[id * 2 + 2];
            }

            @Override
            double getKex(double[] pars) {
                return pars[0];
            }

            @Override
            double getKex(double[] pars, int id) {
                return pars[0];
            }

            @Override
            public int[][] makeMap(int n) {
                int[][] map = new int[n][3];
                for (int i = 0; i < n; i++) {
                    map[i][0] = 0;
                    map[i][1] = 2 * i + 1;
                    map[i][2] = 2 * i + 2;
                }
                return map;
            }

            @Override
            public int[][] makeMap(int n, int m) {
                int[][] map = new int[n][3];
                for (int i = 0; i < n; i++) {
                    map[i][0] = 0;
                    map[i][1] = 2 * i + 1;
                    map[i][2] = 2 * i + 2;
                }
                return map;
            }
        },
        //        ISHIMA("isima", "R2", "Rex", "PaDw", "Tau") {
        //            double calculate(double[] par, double tcp, double field) {
        //                /*
        //                Ishima and Torchia approximation for skewed populations and all time scales
        //                R2(1/tcp)=R2+Rex/(1+Tau*sqrt(144/tcp**4+PaDw**4))
        //                 */
        //                double R2 = par[0];
        //                double Rex = par[1];
        //                double PaDw = par[2];
        //                double Tau = par[3];
        //                double value = R2 + Rex / (1 + Tau * FastMath.sqrt(144.0 / FastMath.pow(tcp, 4) + FastMath.pow(PaDw, 4)));
        //                return value;
        //            }
        //        },
        CPMGSLOW("cpmgslow", 2, "Kex", "pA", "R2", "dW") {
            @Override
            public double calculate(double[] par, int[] map, double nu, int idNum, double field) {
// do we need a different R2 for each field
                double kEx = par[map[0]];
                double pA = par[map[1]]; // p1-p2
                double r2 = par[map[2]];
                double dW = par[map[3]];
                double pB = 1.0 - pA;
                double pDelta = pA - pB;

                double fieldAdjust = field / fieldRef;
                dW *= fieldAdjust;
                dW *= 2.0 * Math.PI;

                double tauCP = 1.0 / (2.0 * nu);
                // fixme should this be plus or minus (cf. Palmer book and NESSY/Guard papers)
                // or next statemetn
                //double psi = kEx * kEx - dW * dW;
                double psi = (pDelta * kEx) * (pDelta * kEx) - dW * dW + 4.0 * pA * pB * kEx * kEx;
                double zeta = -2.0 * dW * kEx * pDelta;
                double eta1 = Math.sqrt(psi * psi + zeta * zeta);
// fixme should the 0.5 be 1.0/sqrt(2.0) ??
                //double etaP = 0.5 * tcp * Math.sqrt(eta1 + psi);
                //double etaM = 0.5 * tcp * Math.sqrt(eta1 - psi);
                double etaP = (1.0 / Math.sqrt(2.0)) * tauCP * Math.sqrt(eta1 + psi);
                double etaM = (1.0 / Math.sqrt(2.0)) * tauCP * Math.sqrt(eta1 - psi);
//System.out.println("pA " + pA + " pB " + pB + " dW " + dW + " kEx " + kEx + " 1/tcp " + (1.0/tcp) + " pd " + pDelta + " zeta " + zeta + " eta " + eta1 + " etaP " + etaP + " etaM " + etaM);
                double d1 = (psi + 2.0 * dW * dW) / Math.sqrt(psi * psi + zeta * zeta);
                // fixme should the sqrt be there (cf. Palmer book and NESSY/Guard papers)
                double dP = 0.5 * (d1 + 1);
                double dM = 0.5 * (d1 - 1);
// should second cos be cosh
                double ch = dP * FastMath.cosh(etaP) - dM * FastMath.cos(etaM);
                double rexContrib = 0.5 * (kEx - (1.0 / tauCP) * FastMath.acosh(ch));
                double dR = (1.0 / tauCP) * FastMath.acosh(ch);
//System.out.println(" d1 " + d1 + " dP " + dP + " dM " + dM + " ch " + ch + " acosh " + FastMath.acosh(ch) + " dR " + dR + " rex " + rexContrib);
                double value = r2 + rexContrib;
                if (Double.isInfinite(value)) {
                    System.out.println("infi");
                    System.out.println(kEx + " pa " + pA + " r2 " + r2 + " dw " + dW + " psi " + psi + " rex " + rexContrib + " ch " + ch + " eta " + etaP + " " + etaM + " " + eta1 + " " + psi);
                }

                return value;

            }

            @Override
            double[] guess(double[] xValues, double[] yValues, int[] idNums, int nID, double field) {
                double[] guesses = new double[2 + nID * 2];
                double kExSum = 0.0;
                double pa = 0.95;
                for (int id = 0; id < nID; id++) {
                    double minY = DataUtil.getMinValue(yValues, idNums, id);
                    double maxY = DataUtil.getMaxValue(yValues, idNums, id);

                    double r2 = minY * 0.8;
                    double rex = maxY - r2;
                    if (rex < 0.0) {
                        rex = 0.0;
                    }
                    double vMid = DataUtil.getMidValue(yValues, xValues, idNums, id);
                    double tauMid = 1.0 / (2.0 * vMid);
                    double kex = 1.915 / (0.5 * tauMid); // 1.915 comes from solving equation iteratively at tcp rex 0.5 half max
                    if (kex > 1000.0) { // fixme what should upper limit be (if kex gets too high can get infinite results
                        kex = 1000.0;
                    }

                    double dw2 = rex / (pa * (1.0 - pa)) * kex;
                    guesses[2 + id * 2] = r2;
                    guesses[3 + id * 2] = Math.sqrt(dw2) / (2.0 * Math.PI);
                    kExSum += kex;
                }
                guesses[0] = kExSum /= nID;
                guesses[1] = pa;
                return guesses;
            }

            @Override
            double[][] boundaries(double[] xValues, double[] yValues, int[] idNums, int nID, double field) {
                double[] guesses = guess(xValues, yValues, idNums, nID, field);
                double[][] boundaries = new double[2][guesses.length];
                boundaries[0][0] = 0.0;
                boundaries[1][0] = guesses[0] * 4;
                boundaries[0][1] = 0.5;
                boundaries[1][1] = 0.99;
                for (int id = 0; id < nID; id++) {
                    boundaries[0][id * 2 + 2] = 0.0;
                    boundaries[1][id * 2 + 2] = guesses[id * 2 + 2] * 4;
                    boundaries[0][id * 2 + 3] = 0.0;
                    boundaries[1][id * 2 + 3] = guesses[id * 2 + 3] * 4;
                }
                return boundaries;
            }
//        CPMGSLOW("cpmgslow", 2, "Kex", "pA", "R2", "dW") {

            @Override
            public double getRex(double[] pars
            ) {
                double rex = 0.0;
                if (pars[3] != 0.0) {
                    // bieri and gooley rex = pa*pb*kex(/(1+(kex/dw)^2))
                    rex = pars[1] * (1.0 - pars[1]) * pars[0] / (1.0 + Math.pow(pars[0] / (2.0 * Math.PI * pars[3]), 2));
                }
                return rex;
            }

            @Override
            public double getRex(double[] pars, int id) {
                double rex = 0.0;
                if (pars[3 + 2 * id] != 0.0) {
                    // bieri and gooley rex = pa*pb*kex(/(1+(kex/dw)^2))
                    rex = pars[1] * (1.0 - pars[1]) * pars[0] / (1.0 + Math.pow(pars[0] / pars[3 + 2 * id], 2));
                }
                return rex;
            }

            @Override
            double getKex(double[] pars
            ) {
                return pars[0];
            }

            @Override
            double getKex(double[] pars, int id) {
                return pars[0];
            }

            @Override
            public int[][] makeMap(int n) {
                int[][] map = new int[n][4];
                for (int i = 0; i < n; i++) {
                    map[i][0] = 0;
                    map[i][1] = 1;
                    map[i][2] = 2 * i + 2;
                    map[i][3] = 2 * i + 3;
                }
                return map;
            }

            @Override
            public int[][] makeMap(int n, int m) {
                int[][] map = new int[n][4];
                for (int i = 0; i < n; i++) {
                    map[i][0] = 0;
                    map[i][1] = 1;
                    map[i][2] = 2 * i + 2;
                    map[i][3] = 2 * i + 3;
                }
                return map;
            }

        }, //        FULLREX("fullrex", "PaPb", "Dw", "kex") {
        //            double calculate(double[] par, double tcp, double field) {
        //                double R2 = par[0];
        //                double Rex = par[1];
        //                double PaDw = par[2];
        //                double Tau = par[3];
        //                double value = 0.0;
        //                return value;
        //            }
        //        },
        //        REX("rex") {
        //            double calculate(double[] par, double tcp, double field) {
        //                /*
        //
        // Fast limit cpmg function without R2 offset
        //
        //     R2(1/tcp)=Rex*(1 - (2*Tau/tcp)*tanh(1/(2*Tau/tcp)))
        //     Function: Rex 
        //     Parameters: Rex, Tau
        //                 */
        //                double R2 = par[0];
        //                double Rex = par[1];
        //                double PaDw = par[2];
        //                double Tau = par[3];
        //                double value = Rex * (1 - (2 * Tau / tcp) * FastMath.tanh(1 / (2 * Tau / tcp)));
        //                return value;
        //            }
        //        },
        //        SITE3("3site") {
        //            double calculate(double[] par, double x, double field) {
        //                /*
        //Fast limit cpmg function for 3 sites
        //     R2(1/tcp)=R2+Sum{Rex*(1 - (2*Tau/tcp)*tanh(1/(2*Tau/tcp)))}
        //     Function: 3-site_CPMG 
        //     Parameters: R2,Rex1, Tau1, Rex2, Tau2
        //                 */
        //                double R2 = par[0];
        //                double Rex1 = par[1];
        //                double Tau1 = par[2];
        //                double Rex2 = par[3];
        //                double Tau2 = par[4];
        //                double tcp = 1.0 / x;
        //
        //                double value = R2;
        //                value += Rex1 * (1 - (2 * Tau1 / tcp) * FastMath.tanh(1 / (2 * Tau1 / tcp)));
        //                value += Rex2 * (1 - (2 * Tau2 / tcp) * FastMath.tanh(1 / (2 * Tau2 / tcp)));
        //                return value;
        //            }
        //        }
        ;

        final String equationName;
        final int nGroupPars;
        String[] parNames;
        double fieldRef;

        public abstract double calculate(double[] par, int[] map, double x, int idNum, double field);

        abstract double[] guess(double[] xValues, double[] yValues, int[] idNums, int nID, double field);

        abstract double[][] boundaries(double[] xValues, double[] yValues, int[] idNums, int nID, double field);

        public abstract double getRex(double[] pars);

        abstract double getRex(double[] pars, int id);

        abstract double getKex(double[] pars);

        abstract double getKex(double[] pars, int id);

        public void setFieldRef(double field) {
            fieldRef = field;
        }

        public abstract int[][] makeMap(int n);

        public abstract int[][] makeMap(int n, int m);

        public String[] getParNames() {
            return parNames;
        }

        public int getNGroupPars() {
            return nGroupPars;
        }

        CPMGEquation(String equationName, int nGroupPars, String... parNames) {
            this.equationName = equationName;
            this.parNames = parNames;
            this.nGroupPars = nGroupPars;
        }

    }

    public CalcRDisp() {
        this.equation = CPMGEquation.CPMGFAST;
    }

    public void setEquation(String eqName) {
        equation = CPMGEquation.valueOf(eqName.toUpperCase());
        /*
        if (eqName.equals("NOEX")) {
            this.equation = CPMGEquation.NOEX;
        } else {
            this.equation = CPMGEquation.CPMGFAST;
        }
         */
    }

    public CalcRDisp(double[] x, double[] y, double[] err, double[] fieldValues, double[] fields) throws IllegalArgumentException {
        this(x, y, err, fieldValues, fields, new int[x.length]);
    }

    public CalcRDisp(double[] x, double[] y, double[] err, double[] fieldValues, double[] fields, int[] idNums) throws IllegalArgumentException {
        this.xValues = x.clone();
        this.yValues = y.clone();
        this.errValues = err.clone();
        this.fieldValues = fieldValues.clone();
        this.fields = fields.clone();
        this.idNums = idNums.clone();
        this.idNums = new int[xValues.length];
        this.equation = CPMGEquation.CPMGFAST;
        this.equation.setFieldRef(fields[0]);
        if (setNID()) {
            throw new IllegalArgumentException("Invalid idNums, some values not used");
        }
    }

    public void setMap(int[][] map) {
        this.map = map;
    }

    public int[][] getMap() {
        return map;
    }

    private boolean setNID() {
        nID = Arrays.stream(idNums).max().getAsInt() + 1;
        boolean[] checkIDs = new boolean[nID];
        Arrays.stream(idNums).forEach(id -> checkIDs[id] = true);
        map = equation.makeMap(nID);
        return IntStream.range(0, checkIDs.length).anyMatch(id -> checkIDs[id] == false);
    }

    public void setXY(double[] x, double[] y) {
        this.xValues = x;
        this.yValues = y;
        this.idNums = new int[x.length];
    }

    public void setErr(double[] err) {
        this.errValues = err;
    }

    public void setFieldValues(double[] fieldValues) {
        this.fieldValues = fieldValues;
    }

    public void setFields(double[] fields) {
        this.fields = fields;
        this.equation.setFieldRef(fields[0]);
    }

    public void setIds(int[] idNums) throws IllegalArgumentException {
        this.idNums = idNums;
        if (setNID()) {
            throw new IllegalArgumentException("Invalid idNums, some values not used");
        }
    }

    public void setAbsMode(boolean value) {
        this.absMode = value;
    }

    public void setNSim(int value) {
        this.nSim = value;
    }

    public class Checker extends SimpleValueChecker {

        public Checker(double relativeThreshold, double absoluteThreshold, int maxIter) {
            super(relativeThreshold, absoluteThreshold, maxIter);
        }

        public boolean converged(final int iteration, final PointValuePair previous, final PointValuePair current) {
            boolean converged = super.converged(iteration, previous, current);
            if (reportFitness) {
                if (converged) {
                    System.out.println(previous.getValue() + " " + current.getValue());
                }
                if (converged || (iteration == 1) || ((iteration % reportAt) == 0)) {
                    long time = System.currentTimeMillis();
                    long deltaTime = time - startTime;
                    System.out.println(deltaTime + " " + iteration + " " + current.getValue());
                }
            }
            return converged;
        }
    }

    // fixme is there a thread safe RandomGenerator
    public final RandomGenerator DEFAULT_RANDOMGENERATOR = new MersenneTwister(1);

    public String[] getParNames() {
        return equation.getParNames();
    }

    public int getNGroupPars() {
        return equation.getNGroupPars();
    }

    @Override
    public double value(double[] par) {
        double sumAbs = 0.0;
        double sumSq = 0.0;
        for (int i = 0; i < xValues.length; i++) {
            final double value;
            value = equation.calculate(par, map[idNums[i]], xValues[i], idNums[i], fieldValues[i]);
            //System.out.println( "xxxxxxxxxxx " + value + " " + yValues[i] + " " + equation.name());
            double delta = (value - yValues[i]);
            //System.out.print(xValues[i] + " " + yValues[i] + " " + value + " " + (delta*delta) + " ");
            //double delta = (value - yValues[i]) / errValues[i];
            sumAbs += FastMath.abs(delta);
            sumSq += delta * delta;
        }
//        if (reportFitness) {
        //           double rms = Math.sqrt(sumSq / xValues.length);
        //          for (double p:par) {
        //             System.out.print(p + " ");
        //        }
        //       System.out.println(" " + sumSq + " " + sumAbs + " " + rms);
        //  }

        if (absMode) {
            return sumAbs;
        } else {
            return sumSq;
        }
    }

    public double getRSS(double[] par) {
        double rss = 0.0;
        for (int i = 0; i < xValues.length; i++) {
            final double value;
            value = equation.calculate(par, map[idNums[i]], xValues[i], idNums[i], fieldValues[i]);
            double delta = value - yValues[i];
            rss += delta * delta;
        }
        return rss;
    }

    public double getRMS(double[] par) {
        double rss = 0.0;
        for (int i = 0; i < xValues.length; i++) {
            final double value;
            value = equation.calculate(par, map[idNums[i]], xValues[i], idNums[i], fieldValues[i]);

            double delta = value - yValues[i];
            rss += delta * delta;
        }
        return Math.sqrt(rss / xValues.length);
    }

    public double getAICc(double[] par) {
        double rss = getRSS(par);
        int k = par.length;
        int n = xValues.length;
        double aic = 2 * k + n * Math.log(rss);
        double aicc = aic + 2 * k * (k + 1) / (n - k - 1);
        return aicc;
    }

    public double[] getPredicted(double[] par) {
        double[] yPred = new double[yValues.length];
        double rss = 0.0;
        for (int i = 0; i < xValues.length; i++) {
            yPred[i] = equation.calculate(par, map[idNums[i]], xValues[i], idNums[i], fieldValues[i]);
        }
        return yPred;
    }

    public void dump(double[] par) {
        for (double field : fields) {
            System.out.println("# field " + field);
            for (int i = 0; i < xValues.length; i++) {
                if (FastMath.abs(fieldValues[i] - field) < 0.01) {
                    double yCalc = equation.calculate(par, map[idNums[i]], xValues[i], idNums[i], field);
                    System.out.printf("%8.5f %8.5f %8.5f\n", xValues[i], yValues[i], yCalc);
                }
            }
        }

        DescriptiveStatistics dstat = new DescriptiveStatistics(xValues);
        double min = dstat.getMin();
        double max = dstat.getMax();
        int n = 20;
        double delta = (max - min) / (n - 1);
        for (double field : fields) {
            System.out.println("# field " + field);
            for (int i = 0; i < n; i++) {
                double x = min + i * delta;
                double y = equation.calculate(par, map[idNums[i]], x, idNums[i], field);
                System.out.printf("%8.5f %8.5f\n", x, y);
            }
        }

    }

    public void dump(double[] par, double[] xValues, double field) {
        equation.setFieldRef(field);
        for (int i = 0; i < xValues.length; i++) {
            double yCalc = equation.calculate(par, map[idNums[i]], xValues[i], idNums[i], field);
            System.out.printf("%8.5f %8.5f %8.5f\n", xValues[i], yCalc, field);
        }
    }

    public ArrayList<Double> simY(double[] par, double[] xValues, double field) {
        ArrayList<Double> result = new ArrayList<>();
        equation.setFieldRef(field);
        for (int i = 0; i < xValues.length; i++) {
            double yCalc = equation.calculate(par, map[idNums[i]], xValues[i], idNums[i], field);
            result.add(yCalc);
        }
        return result;
    }

    public static double rDisp(double[] par, double x) {
        double a = par[0];
        double b = par[1];
        double c = par[2];
        // fixme check for x=0;
        return Math.exp(-(a - a * Math.sin(b / x) / (b / x) + c));
    }

    public double[] guess() {
        double[] guess = equation.guess(xValues, yValues, idNums, nID, fieldValues[0]);
        return guess;
    }

    public double[][] boundaries() {
        //return equation.boundaries(xValues, yValues, fieldValues[0]);
        double[][] boundaries = equation.boundaries(xValues, yValues, idNums, nID, fieldValues[0]);
        return boundaries;
    }

    public double[] getRex(double[] pars) {
        double[] result = new double[nID];
        for (int i = 0; i < nID; i++) {
            result[i] = equation.getRex(pars, i);
        }
        return result;
    }

    public double getRex(double[] pars, int id) {
        return equation.getRex(pars, id);
    }

    public double[] getRexError() {
        return rexErrors.clone();
    }

    public double getKex(double[] pars) {
        return equation.getKex(pars);
    }

    public PointValuePair refine(double[] guess, double[] lowerBounds, double[] upperBounds, double[] inputSigma) {
        startTime = System.currentTimeMillis();
        DEFAULT_RANDOMGENERATOR.setSeed(1);
        double lambdaMul = 3.0;
        int lambda = (int) (lambdaMul * FastMath.round(4 + 3 * FastMath.log(guess.length)));
        //int nSteps = guess.length*1000;
        int nSteps = 2000;
        double stopFitness = 0.0;
        int diagOnly = 0;
        //new Checker(100 * Precision.EPSILON, 100 * Precision.SAFE_MIN, nSteps));
        CMAESOptimizer optimizer = new CMAESOptimizer(nSteps, stopFitness, true, diagOnly, 0,
                DEFAULT_RANDOMGENERATOR, true,
                new Checker(1.0e-7, 1.0e-7, nSteps));
        PointValuePair result = null;

        try {
            result = optimizer.optimize(
                    new CMAESOptimizer.PopulationSize(lambda),
                    new CMAESOptimizer.Sigma(inputSigma),
                    new MaxEval(2000000),
                    new ObjectiveFunction(this), GoalType.MINIMIZE,
                    new SimpleBounds(lowerBounds, upperBounds),
                    new InitialGuess(guess));
        } catch (Exception e) {
            e.printStackTrace();
        }
        return result;

    }

    public double[] simBounds(double[] start, double[] lowerBounds, double[] upperBounds, double[] inputSigma) {
        reportFitness = false;
        int nPar = start.length;
        double[][] parValues = new double[nPar][nSim];
        double[] yPred = getPredicted(start);
        double[] yValuesOrig = yValues.clone();
        double[][] rexValues = new double[nID][nSim];
        rexErrors = new double[nID];
        for (int i = 0; i < nSim; i++) {
            for (int k = 0; k < yValues.length; k++) {
                yValues[k] = yPred[k] + errValues[k] * random.nextGaussian();
            }
            PointValuePair result = refine(start, lowerBounds, upperBounds, inputSigma);
            double[] rPoint = result.getPoint();
            for (int j = 0; j < nPar; j++) {
                parValues[j][i] = rPoint[j];
            }
            if (equation == CPMGEquation.CPMGSLOW) {
                for (int j = 0; j < nID; j++) {
                    rexValues[j][i] = equation.getRex(result.getPoint(), j);
                }
            }
        }
        double[] parSDev = new double[nPar];
        for (int i = 0; i < nPar; i++) {
            DescriptiveStatistics dStat = new DescriptiveStatistics(parValues[i]);
            double p5 = dStat.getPercentile(5.0);
            double p95 = dStat.getPercentile(95.0);
            parSDev[i] = dStat.getStandardDeviation();
        }
        if (equation == CPMGEquation.CPMGSLOW) {
            for (int j = 0; j < nID; j++) {
                DescriptiveStatistics dStat = new DescriptiveStatistics(rexValues[j]);
                rexErrors[j] = dStat.getStandardDeviation();
            }
        }
        yValues = yValuesOrig;
        return parSDev;
    }

    public double[] simBoundsStream(double[] start, double[] lowerBounds, double[] upperBounds, double[] inputSigma) {
        reportFitness = false;
        int nPar = start.length;
        double[][] parValues = new double[nPar][nSim];
        double[][] rexValues = new double[nID][nSim];
        rexErrors = new double[nID];
        double[] yPred = getPredicted(start);
        IntStream.range(0, nSim).parallel().forEach(i -> {
            CalcRDisp rDisp = new CalcRDisp(xValues, yPred, errValues, fieldValues, fields, idNums);
            rDisp.setEquation(equation.name());
            rDisp.setAbsMode(absMode);
            double[] newY = new double[yValues.length];
            for (int k = 0; k < yValues.length; k++) {
                newY[k] = yPred[k] + errValues[k] * random.nextGaussian();
            }
            rDisp.setXY(xValues, newY);
            rDisp.setIds(idNums);

            PointValuePair result = rDisp.refine(start, lowerBounds, upperBounds, inputSigma);
            double[] rPoint = result.getPoint();
            for (int j = 0; j < nPar; j++) {
                parValues[j][i] = rPoint[j];
            }
            if (equation == CPMGEquation.CPMGSLOW) {
                for (int j = 0; j < nID; j++) {
                    rexValues[j][i] = equation.getRex(result.getPoint(), j);
                }
            }
        });

        double[] parSDev = new double[nPar];
        for (int i = 0; i < nPar; i++) {
            DescriptiveStatistics dStat = new DescriptiveStatistics(parValues[i]);
            parSDev[i] = dStat.getStandardDeviation();
        }
        if (equation == CPMGEquation.CPMGSLOW) {
            for (int j = 0; j < nID; j++) {
                DescriptiveStatistics dStat = new DescriptiveStatistics(rexValues[j]);
                rexErrors[j] = dStat.getStandardDeviation();
            }
        }
        return parSDev;
    }

    public double[] simBoundsBootstrapStream(double[] start, double[] lowerBounds, double[] upperBounds, double[] inputSigma) {
        reportFitness = false;
        int nPar = start.length;
        double[][] parValues = new double[nPar][nSim];
        double[][] rexValues = new double[nID][nSim];
        rexErrors = new double[nID];
        IntStream.range(0, nSim).parallel().forEach(i -> {
            CalcRDisp rDisp = new CalcRDisp(xValues, yValues, errValues, fieldValues, fields, idNums);
            rDisp.setEquation(equation.name());
            rDisp.setAbsMode(absMode);
            double[] newX = new double[yValues.length];
            double[] newY = new double[yValues.length];
            double[] newErr = new double[yValues.length];
            double[] newFieldValues = new double[yValues.length];
            int[] newID = new int[yValues.length];
            for (int k = 0; k < yValues.length; k++) {
                int rI = random.nextInt(yValues.length);
                newX[k] = xValues[rI];
                newY[k] = yValues[rI];
                newErr[k] = errValues[rI];
                newFieldValues[k] = fieldValues[rI];
                newID[k] = idNums[rI];
            }
            // fixme  idNum should be set in above loop
            rDisp.setXY(newX, newY);
            rDisp.setErr(newErr);
            rDisp.setFieldValues(newFieldValues);
            rDisp.setIds(newID);

            PointValuePair result = rDisp.refine(start, lowerBounds, upperBounds, inputSigma);
            double[] rPoint = result.getPoint();
            for (int j = 0; j < nPar; j++) {
                parValues[j][i] = rPoint[j];
            }
            if (equation == CPMGEquation.CPMGSLOW) {
                for (int j = 0; j < nID; j++) {
                    rexValues[j][i] = equation.getRex(result.getPoint(), j);
                }
            }
        });

        double[] parSDev = new double[nPar];
        for (int i = 0; i < nPar; i++) {
            DescriptiveStatistics dStat = new DescriptiveStatistics(parValues[i]);
            parSDev[i] = dStat.getStandardDeviation();
        }
        if (equation == CPMGEquation.CPMGSLOW) {
            for (int j = 0; j < nID; j++) {
                DescriptiveStatistics dStat = new DescriptiveStatistics(rexValues[j]);
                rexErrors[j] = dStat.getStandardDeviation();
            }
        }
        return parSDev;
    }

}

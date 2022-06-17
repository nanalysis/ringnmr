package org.comdnmr.modelfree;

import java.util.ArrayList;
import java.util.List;

public class SpectralDensityCalculator {

    public static double[][] calcJR1R2NOE(List<RelaxDataValue> dataValues) {
        int nDataValues = dataValues.size();
        int nFreq = 1 + 2 * nDataValues;
        double[][] result = new double[3][nFreq];
        int iField = 0;
        for (RelaxDataValue value : dataValues) {
            R1R2NOEDataValue relaxDataValue = (R1R2NOEDataValue) value;
            RelaxEquations relaxEq = relaxDataValue.relaxObj;
            double r1 = relaxDataValue.R1;
            double r2 = relaxDataValue.R2;
            double noe = relaxDataValue.NOE;
            double r1Err = relaxDataValue.R1err;
            double r2Err = relaxDataValue.R2err;
            double noeErr = relaxDataValue.NOEerr;
            double sigma = (noe - 1.0) * r1 * RelaxEquations.GAMMA_N / RelaxEquations.GAMMA_H;
            double sigmaErr = sigma * Math.sqrt(Math.pow((noeErr / (noe - 1.0)), 2) + Math.pow((r1Err / r1), 2));

            double d2 = relaxEq.getD2();
            double c2 = relaxEq.getC2();

            double j87H = 4.0 * sigma / (5.0 * d2);
            double j87Herr = 4.0 * sigmaErr / (5.0 * d2);

            double jNMul = 4.0 / (3.0 * d2 + 4.0 * c2);
            double jN = (r1 - 1.249 * sigma) * jNMul;
            double jNerr = jNMul * Math.sqrt(Math.pow(r1Err, 2) + Math.pow(1.249 * sigmaErr, 2));

            double j0Mul = 6.0 / (3.0 * d2 + 4.0 * c2);
            double j0 = j0Mul * (r2 - 0.5 * r1 - 0.454 * sigma);
            double j0Err = j0Mul * Math.sqrt(Math.pow(r2Err, 2) + Math.pow(0.5 * r1Err, 2) + Math.pow(0.454 * sigmaErr, 2));

            result[0][0] = 0.0;
            result[1][0] += j0;
            result[2][0] += j0Err * j0Err;

            result[0][iField * 2 + 1] = 0.87 * relaxEq.getWI();
            result[1][iField * 2 + 1] = j87H;
            result[2][iField * 2 + 1] = j87Herr;

            result[0][iField * 2 + 2] = relaxEq.getWS();
            result[1][iField * 2 + 2] = jN;
            result[2][iField * 2 + 2] = jNerr;
            iField++;
        }
        result[1][0] /= nDataValues;
        result[2][0] = Math.sqrt(result[2][0]);
        return result;
    }

    public static double[][] calcJDeuterium(List<RelaxDataValue> dataValues) {
        double[][] result = null;
        if (!dataValues.isEmpty()) {
            List<Double> rValues = new ArrayList<>();
            List<Double> errValues = new ArrayList<>();
            List<Double> fields = new ArrayList<>();
            for (var value : dataValues) {
                var dValue = (DeuteriumDataValue) value;
                rValues.add(dValue.R1);
                rValues.add(dValue.R2);
                rValues.add(dValue.rQ);
                rValues.add(dValue.rAP);
                errValues.add(dValue.R1err);
                errValues.add(dValue.R2err);
                errValues.add(dValue.rQError);
                errValues.add(dValue.rAPError);
                fields.add(dValue.relaxObj.getSF() * RelaxEquations.GAMMA_D / RelaxEquations.GAMMA_H * 2.0 * Math.PI);
            }
            result = DeuteriumMapping.jointMapping(rValues, errValues, fields);
        }
        return result;
    }
}

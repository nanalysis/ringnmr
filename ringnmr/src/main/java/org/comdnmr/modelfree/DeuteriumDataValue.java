package org.comdnmr.modelfree;

import org.comdnmr.modelfree.MolDataValues;
import org.comdnmr.modelfree.RelaxDataValue;
import org.comdnmr.modelfree.RelaxEquations;

public class DeuteriumDataValue extends RelaxDataValue {
    final double rQ;
    final double rQError;
    final double rAP;
    final double rAPError;

    public DeuteriumDataValue(MolDataValues molDataValue, double r1,
                              double r1Error, double r2, double r2Error,
                              double rQ, double rQError, double rAP, double rAPError,
                              RelaxEquations relaxObj) {
        super(molDataValue, r1, r1Error, r2, r2Error, relaxObj);
        this.rAP = rAP;
        this.rAPError = rAPError;
        this.rQ = rQ;
        this.rQError = rQError;
    }
}

package org.comdnmr.modelfree;

import org.apache.commons.math3.optim.PointValuePair;
import org.comdnmr.modelfree.models.MFModelIso;

import java.util.Map;
import java.util.Random;

public class FitModel {
    Double tau;
    boolean fitTau = false;
    boolean fitJ = false;
    boolean fitExchange = false;
    double tauFraction = 0.25;
    double lambda = 0.0;
    double t2Limit = 0.0;
    int nReplicates = 0;

    double[][] replicates(Map<String, MolDataValues> molDataRes,
                          MFModelIso bestModel, double localTauFraction,
                          boolean localFitTau, double[] pars, Random random) {
        double[][] repData = new double[pars.length][nReplicates];
        for (int iRep = 0; iRep < nReplicates; iRep++) {
            Score score2 = fitReplicate(molDataRes, bestModel, localTauFraction, localFitTau, pars, random);
            double[] repPars = score2.getPars();
            for (int iPar = 0; iPar < pars.length; iPar++) {
                repData[iPar][iRep] = repPars[iPar];
            }
        }
        return repData;
    }

    Score fitReplicate(Map<String, MolDataValues> molDataRes, MFModelIso model,
                       double localTauFraction, boolean localFitTau, double[] pars, Random random) {
        RelaxFit relaxFit = new RelaxFit();
        relaxFit.setRelaxData(molDataRes);
        relaxFit.setLambda(lambda);
        relaxFit.setFitJ(fitJ);
        Map<String, MolDataValues> molDataMap = relaxFit.genBootstrap(random, model, pars);
        relaxFit.setRelaxData(molDataMap);

        model.setTauFraction(localTauFraction);
        double[] lower = model.getLower();
        double[] upper = model.getUpper();
        PointValuePair fitResult = relaxFit.fitResidueToModel(pars, lower, upper);
        return relaxFit.score(fitResult.getPoint(), true);
    }

    public Double getTau() {
        return tau;
    }

    public void setTau(Double value) {
        tau = value;
    }

    public void setFitTau(boolean value) {
        fitTau = value;
    }

    public void setLambda(double value) {
        this.lambda = value;
    }

    public void setFitJ(boolean value) {
        this.fitJ = value;
    }

    public void setNReplicates(int value) {
        this.nReplicates = value;
    }

    public void setT2Limit(double value) {
        this.t2Limit = value;
    }
}

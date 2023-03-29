package org.comdnmr.modelfree;

import javafx.beans.property.ReadOnlyObjectProperty;
import javafx.concurrent.Service;
import javafx.concurrent.Task;
import javafx.concurrent.Worker;
import org.apache.commons.math3.optim.PointValuePair;
import org.comdnmr.data.ExperimentSet;
import org.comdnmr.modelfree.models.MFModelIso;
import org.comdnmr.util.CoMDOptions;
import org.comdnmr.util.ProcessingStatus;
import org.nmrfx.chemistry.relax.ResonanceSource;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.function.Function;

public abstract class FitModel {
    Double tau;
    boolean fitTau = false;
    boolean fitJ = false;
    boolean bootstrap = false;
    boolean fitExchange = false;
    double tauFraction = 0.25;
    double lambda = 0.0;
    boolean useLambda = false;
    double t2Limit = 0.0;
    int nReplicates = 0;
    private final FitResidues fitResidues = new FitResidues();
    final ReadOnlyObjectProperty<Worker.State> stateProperty = fitResidues.worker.stateProperty();
    private boolean isProcessing = false;
    List<String> modelNames = new ArrayList<>();
    String searchKey = null;
    AtomicBoolean cancelled = new AtomicBoolean(false);

    Function<Double, Double> updaterFunction;
    Function<ProcessingStatus, Double> statusFunction;


    public void setup(String searchKey, List<String> modelNames) {
        this.searchKey = searchKey;
        this.modelNames.clear();
        this.modelNames.addAll(modelNames);
    }

    public void updaters(Function<Double, Double> updaterFunction, Function<ProcessingStatus, Double> statusFunction) {
        this.updaterFunction = updaterFunction;
        this.statusFunction = statusFunction;
    }

    public abstract void testIsoModel();

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
        relaxFit.setUseLambda(useLambda);
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

    public void setUseLambda(boolean value) {
        this.useLambda = value;
    }

    public boolean useLambda() {
        return useLambda && lambda > 1.0e-8;
    }

    public void setFitJ(boolean value) {
        this.fitJ = value;
    }

    public void setBootstrap(boolean value) {
        this.bootstrap = value;
    }

    public void setNReplicates(int value) {
        this.nReplicates = value;
    }

    public void setT2Limit(double value) {
        this.t2Limit = value;
    }

    public void setTauFraction(double value) {
        tauFraction = value;
    }

    private class FitResidues {

        String script;
        public Worker<Integer> worker;

        private FitResidues() {
            worker = new Service<Integer>() {

                @Override
                protected Task createTask() {
                    return new Task() {
                        @Override
                        protected Object call() {
                            cancelled.set(false);
                            updateStatus("Start processing");
                            updateTitle("Start Processing");
                            fitAll(this);
                            return 0;
                        }
                    };
                }
            };

            ((Service<Integer>) worker).setOnSucceeded(event -> {
                finishProcessing();
                setProcessingOff();
            });
            ((Service<Integer>) worker).setOnCancelled(event -> {
                cancelled.set(true);
                setProcessingOff();
                setProcessingStatus("cancelled", false);
            });
            ((Service<Integer>) worker).setOnFailed(event -> {
                setProcessingOff();
                final Throwable exception = worker.getException();
                setProcessingStatus(exception.getMessage(), false, exception);

            });

        }
    }


    public void fitResidues() {
        setProcessingOn();
        updateProgress(0.0);
        ((Service) fitResidues.worker).restart();
    }

    public void updateProgress(Double f) {
        if (updaterFunction != null) {
            updaterFunction.apply(f);
        }
    }

    public void updateStatus(String s) {
        setProcessingStatus(s, true);
    }

    public void setProcessingStatus(String s, boolean ok) {
        setProcessingStatus(s, ok, null);
    }

    public void setProcessingStatus(String s, boolean ok, Throwable throwable) {
        if (statusFunction != null) {
            ProcessingStatus status = new ProcessingStatus(s, ok, throwable);
            statusFunction.apply(status);
        }
    }

    public void clearProcessingTextLabel() {
        ProcessingStatus status = new ProcessingStatus("");
        statusFunction.apply(status);
    }

    public void haltFit() {
        fitResidues.worker.cancel();
    }

    synchronized void setProcessingOn() {
        isProcessing = true;
    }

    synchronized void setProcessingOff() {
        isProcessing = false;
    }

    boolean isProcessing() {
        return isProcessing;
    }

    void finishProcessing() {
        updateStatus("Done");
        try {
        } catch (IllegalArgumentException iAE) {
        }
    }

    public void fitAll(Task task) {
        testIsoModel();
        int nFit = 0;
        int nGroups = 1;
        nFit++;
        if (task != null) {
            updateProgress((1.0 * nFit) / nGroups);
        }
    }
}


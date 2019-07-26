package org.comdnmr.fit.calc;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.function.Function;
import javafx.beans.property.ReadOnlyObjectProperty;
import javafx.concurrent.Service;
import javafx.concurrent.Task;
import javafx.concurrent.Worker;

/**
 *
 * @author Bruce Johnson
 */
public class ResidueFitter {

    private final FitResidues processDataset = new FitResidues();
    final ReadOnlyObjectProperty<Worker.State> stateProperty = processDataset.worker.stateProperty();
    private boolean isProcessing = false;
    ResidueProperties resProps;
    Function<Double, Double> updaterFunction;
    Function<ProcessingStatus, Double> statusFunction;
    List<List<String>> residueFitGroups = null;
    CPMGFitResult fitResult;

    public ResidueFitter() {
    }

    public ResidueFitter(Function<Double, Double> updaterFunction, Function<ProcessingStatus, Double> statusFunction) {
        this.updaterFunction = updaterFunction;
        this.statusFunction = statusFunction;
    }

    public void fitResidues(ResidueProperties resProps) {
        fitResidues(resProps, null);
    }

    public void fitResidues(ResidueProperties resProps, List<List<String>> residueFitGroups) {
        this.resProps = resProps;
        resProps.setupMaps();
        this.residueFitGroups = residueFitGroups;
        setProcessingOn();
        updateProgress(0.0);
        if (residueFitGroups == null) {
            resProps.clearResidueMap();
        }
        ((Service) processDataset.worker).restart();
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
        processDataset.worker.cancel();
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
    }

    public void fitAllResidueGroups(Task task) {
        for (List<String> resGroup : residueFitGroups) {
            int nResidues = resGroup.size();
            int nFit = 0;
            if (task.isCancelled()) {
                break;
            }
            String[] resNumGroup = new String[resGroup.size()];
            resGroup.toArray(resNumGroup);
            List<ResidueInfo> resInfoList = fitResidues(resProps, resNumGroup, nFit, null);
            for (ResidueInfo resInfo : resInfoList) {
                int fitResNum = resInfo.resNum;
                resProps.addResidueInfo(String.valueOf(fitResNum), resInfo);
            }
            nFit++;
            updateProgress((1.0 * nFit) / nResidues);
        }
    }

    public List<List<String>> getAllResidues() {
        Map<String, ExperimentData> expDataSets = resProps.expMaps;
        Set<Integer> resNums = new TreeSet<>();
        for (ExperimentData expData : expDataSets.values()) {
            for (String resNumS : expData.getResidues()) {
                resNums.add(Integer.parseInt(resNumS));
            }
        }
        List<List<String>> allResidues = new ArrayList<>();
        for (Integer resNum : resNums) {
            int groupIndex = -1;
            String resStr = String.valueOf(resNum);
            if (residueFitGroups != null) {
                int i = 0;
                for (List<String> resGroup : residueFitGroups) {
                    if (resGroup.contains(resStr)) {
                        groupIndex = i;
                        break;
                    }
                    i++;
                }
            }
            if (groupIndex != -1) {
                List<String> resGroup = residueFitGroups.get(groupIndex);
                if (resGroup.get(0).equals(resStr)) {
                    allResidues.add(resGroup);
                }
            } else if (residueFitGroups == null) {
                List<String> singleRes = new ArrayList<>();
                singleRes.add(resStr);
                allResidues.add(singleRes);
            }
        }
        return allResidues;
    }

    public void fitAllResidues(Task task) {
        List<List<String>> allResidues = getAllResidues();
        int nGroups = allResidues.size();
        int nFit = 0;
        for (List<String> resList : allResidues) {
            if (task.isCancelled()) {
                break;
            }
            String[] resNumGroup = new String[resList.size()];
            resList.toArray(resNumGroup);
            List<ResidueInfo> resInfoList = fitResidues(resProps, resNumGroup, nFit, null);
            for (ResidueInfo resInfo : resInfoList) {
                int fitResNum = resInfo.resNum;
                resProps.addResidueInfo(String.valueOf(fitResNum), resInfo);
            }
            nFit++;
            updateProgress((1.0 * nFit) / nGroups);

        }

    }

    EquationFitter getFitter() {
        EquationFitter fitter;
        switch (resProps.getExpMode()) {
            case "cpmg":
                fitter = new CPMGFit();
                break;
            case "exp":
                fitter = new ExpFit();
                break;
            case "cest":
                fitter = new CESTFit();
                break;
            case "r1rho":
                fitter = new R1RhoFit();
                break;
            default:
                throw new IllegalArgumentException("Invalid mode " + resProps.getExpMode());
        }
        return fitter;
    }

    public List<ResidueInfo> doNOE(ResidueProperties resProps, String[] resNums, int groupId, String useEquation) {
        List<ResidueInfo> resInfoList = new ArrayList<>();
        Map<String, ResidueInfo> resMap = new HashMap<>();
        for (String residueNumber : resNums) {
            ResidueInfo residueInfo = new ResidueInfo(resProps, Integer.parseInt(residueNumber), groupId, resNums.length);
            resInfoList.add(residueInfo);
            resMap.put(residueNumber, residueInfo);
            residueInfo.addFitResult(null);
        }
        return resInfoList;

    }

    public List<ResidueInfo> fitResidues(ResidueProperties resProps, String[] resNums, int groupId, String useEquation) {
        this.resProps = resProps;
        resProps.setupMaps();
        Map<String, CPMGFitResult> fitResults = new HashMap<>();
        double aicMin = Double.MAX_VALUE;
        String bestEquation = "NOEX";
        List<String> equationNames;
        switch (resProps.getExpMode()) {
            case "cpmg":
                equationNames = CPMGFit.getEquationNames();
                break;
            case "exp":
                equationNames = ExpFit.getEquationNames();
                bestEquation = "EXPAB";
                break;
            case "cest":
                equationNames = CESTFit.getEquationNames();
                bestEquation = "CESTR1RHOPERTURBATIONNOEX";
                break;
            case "r1rho":
                equationNames = R1RhoFit.getEquationNames();
                bestEquation = "R1RHOPERTURBATIONNOEX";
                break;
            case "noe":
                return doNOE(resProps, resNums, groupId, useEquation);
            default:
                throw new IllegalArgumentException("Invalid mode " + resProps.getExpMode());
        }
        for (String equationName : equationNames) {
            if ((useEquation != null) && !equationName.equals(useEquation)) {
                continue;
            }

            EquationFitter equationFitter = getFitter();
            equationFitter.setData(resProps, resNums);
            fitResult = equationFitter.doFit(equationName, null);
            fitResults.put(equationName, fitResult);
            if (fitResult.getAicc() < aicMin) {
                aicMin = fitResult.getAicc();
                if (fitResult.exchangeValid()) {
                    bestEquation = equationName;
                }
            }
//            System.out.println("fit " + fitResult.getAicc() + " " + aicMin + " " + bestEquation);
        }
        List<ResidueInfo> resInfoList = new ArrayList<>();
        Map<String, ResidueInfo> resMap = new HashMap<>();
        for (String residueNumber : resNums) {
            ResidueInfo residueInfo = new ResidueInfo(resProps, Integer.parseInt(residueNumber), groupId, resNums.length);
            resInfoList.add(residueInfo);
            resMap.put(residueNumber, residueInfo);
            for (String equationName : equationNames) {
                residueInfo.addFitResult(fitResults.get(equationName));
            }
        }

        for (String equationName : equationNames) {
            if ((useEquation != null) && !equationName.equals(useEquation)) {
                continue;
            }
            fitResult = fitResults.get(equationName);

            int nCurves = fitResult.getNCurves();
            for (int iCurve = 0; iCurve < nCurves; iCurve++) {
                CurveFit curveFit = fitResult.getCurveFit(iCurve);
                String residueNumber = curveFit.getResNum();
                ResidueInfo residueInfo = resMap.get(residueNumber);
                residueInfo.addCurveSet(curveFit, bestEquation.equals(equationName));
            }
        }
        return resInfoList;
    }

    public CPMGFitResult getFitResult() {
        return fitResult;
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
                            updateStatus("Start processing");
                            updateTitle("Start Processing");
                            fitAllResidues(this);
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

    public static EquationType getEquationType(String name) throws IllegalArgumentException {
        EquationType equationType;

        try {
            equationType = CPMGEquation.valueOf(name);
        } catch (IllegalArgumentException iaE) {
            try {
                equationType = ExpEquation.valueOf(name);
            } catch (IllegalArgumentException iaE2) {
                try {
                    equationType = CESTEquation.valueOf(name);
                } catch (IllegalArgumentException iaE3) {
                    equationType = R1RhoEquation.valueOf(name);
                }
            }
        }
        return equationType;
    }

}

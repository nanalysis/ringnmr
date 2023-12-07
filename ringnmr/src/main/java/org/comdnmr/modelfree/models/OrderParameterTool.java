package org.comdnmr.modelfree.models;

import org.comdnmr.modelfree.FitDeuteriumModel;
import org.comdnmr.modelfree.MolDataValues;
import org.nmrfx.chemistry.Atom;
import org.nmrfx.chemistry.relax.SpectralDensity;

import java.util.Comparator;
import java.util.Map;

public class OrderParameterTool {
    public static void calculateOrderParameters() {
        var data = FitDeuteriumModel.getData(false);
        data.entrySet().stream().sorted(Comparator.comparing(Map.Entry::getKey)).forEach(e -> {
            MolDataValues resData = e.getValue();
            Atom atom = resData.getAtom();
            String key = e.getKey();
            double[][] jValues = resData.getJValues();
            SpectralDensity spectralDensity = new SpectralDensity(key, jValues);
            atom.addSpectralDensity(key, spectralDensity);
        });
    }
}

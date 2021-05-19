/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.modelfree;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.optim.PointValuePair;
import org.comdnmr.modelfree.models.MFModelIso;
import org.comdnmr.modelfree.models.MFModelIso1;
import org.nmrfx.chemistry.Atom;
import org.nmrfx.chemistry.MoleculeBase;
import org.nmrfx.chemistry.MoleculeFactory;
import org.nmrfx.chemistry.Polymer;
import org.nmrfx.chemistry.RelaxationData;
import org.nmrfx.chemistry.RelaxationRex;
import org.nmrfx.chemistry.Residue;

/**
 *
 * @author brucejohnson
 */
public class FitModel {

    Map<String, MolDataValues> getData() {
        Map<String, MolDataValues> molDataValues = new HashMap<>();
        MoleculeBase mol = MoleculeFactory.getActive();
        if (mol != null) {
            for (Polymer polymer : mol.getPolymers()) {
                for (Residue residue : polymer.getResidues()) {
                    var hAtom = residue.getAtom("H");
                    var nAtom = residue.getAtom("N");
                    if ((hAtom != null) && (nAtom != null)) {
                        var hPoint = hAtom.getPoint();
                        var nPoint = nAtom.getPoint();
                        if ((hPoint != null) && (nPoint != null)) {
                            var relaxData = hAtom.getRelaxationData();
                            if ((relaxData == null) || relaxData.isEmpty()) {
                                relaxData = nAtom.getRelaxationData();
                            }
                            if (relaxData != null) {
                                Vector3D vec = hPoint.subtract(nPoint);
                                MolDataValues molData = new MolDataValues(nAtom, vec.toArray());
                                Set<Long> fields = new HashSet<>();
                                for (var entry : relaxData.entrySet()) {
                                    RelaxationData data = entry.getValue();
                                    fields.add(Math.round(data.getField()));
                                }
                                for (var field : fields) {
                                    Double r1 = null;
                                    Double r1Error = null;

                                    Double r2 = null;
                                    Double r2Error = null;

                                    Double noe = null;
                                    Double noeError = null;

                                    for (var entry : relaxData.entrySet()) {
                                        RelaxationData data = entry.getValue();
                                        if (Math.round(data.getField()) == field) {
                                            switch (data.getExpType()) {
                                                case T1:
                                                    r1 = data.getValue();
                                                    r1Error = data.getError();
                                                    break;
                                                case T2:
                                                    r2 = data.getValue();
                                                    r2Error = data.getError();
                                                    break;
                                                case NOE:
                                                    noe = data.getValue();
                                                    noeError = data.getError();
                                                    break;
                                            }
                                        }
                                    }
                                    if ((r1 != null) && (r2 != null) && (noe != null)) {
                                        String e1 = hAtom.getElementName();
                                        String e2 = nAtom.getElementName();
                                        RelaxEquations relaxObj = RelaxEquations.getRelaxEquations(field * 1e6, "H", "N");

                                        RelaxDataValue dValue = new RelaxDataValue(molData, r1, r1Error, r2, r2Error, noe, noeError, relaxObj);
                                        molData.addData(dValue);
                                    }
                                }
                                if (!molData.dataValues.isEmpty()) {
                                    molDataValues.put(molData.specifier, molData);
                                }
                            }
                        }
                    }

                }
            }
        }
        return molDataValues;
    }

    public void testIsoModel() {
        Map<String, MolDataValues> molData = getData();
        if (!molData.isEmpty()) {
            double tau = 37.5e-9;
            MFModelIso model = new MFModelIso1(tau, true);
            testIsoModel(model, molData, null);
        }
    }

    public void testIsoModel(MFModelIso model, Map<String, MolDataValues> molData, String matchSpec) {
        RelaxFit relaxFit = new RelaxFit();
        Map<String, MolDataValues> molDataRes = new TreeMap<>();
        double tau = 36.0e-9;
        Map<String, String> extras = new HashMap<>();
        MoleculeBase mol = MoleculeFactory.getActive();

        for (String key : molData.keySet()) {
            if ((matchSpec != null) && !key.equals(matchSpec)) {
                continue;
            }
            molDataRes.clear();
            MolDataValues resData = molData.get(key);
            if (!resData.getData().isEmpty()) {
                resData.setTestModel(model);
                molDataRes.put(key, molData.get(key));
                relaxFit.setRelaxData(molDataRes);
                double[] start = model.getStart(tau, false);
                double[] lower = model.getLower(tau, false);
                double[] upper = model.getUpper(tau, false);
                PointValuePair fitResult = relaxFit.fitResidueToModel(start, lower, upper);
                double[] values = fitResult.getPoint();
                double score = fitResult.getValue();
                double orderPar = values[0];
                double rexValue = values[1];
                RelaxationRex relaxData = new RelaxationRex("order", RelaxationData.relaxTypes.S2,
                        new ArrayList<>(), 700.0, 298, orderPar, 0.0,
                        rexValue, 0.0, extras);
                Atom atom = mol.findAtom(key);
                atom.relaxData.put("order", relaxData);

            }
        }
    }

}

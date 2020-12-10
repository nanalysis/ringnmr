/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.data;

import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.nmrfx.chemistry.Atom;
import org.nmrfx.chemistry.MoleculeBase;
import org.nmrfx.chemistry.MoleculeFactory;
import org.nmrfx.chemistry.Polymer;
import org.nmrfx.chemistry.Residue;
import org.nmrfx.chemistry.io.PDBAtomParser;
import org.nmrfx.peaks.Peak;
import org.nmrfx.peaks.PeakList;
import org.nmrfx.peaks.SpectralDim;

/**
 *
 * @author brucejohnson
 */
public class DynamicsSource {

    private static final String RESIDUE_STR = "^(([a-zA-Z]+):)?(([a-zA-Z])?(-?[0-9]+)|([ij]))(\\.([a-zA-Z].*))";
    private static final String PEAK_STR = "^([a-zA-Z].*)\\.([0-9]+)";
    private static final Pattern RESIDUE_PATTERN = Pattern.compile(RESIDUE_STR);
    private static final Pattern PEAK_PATTERN = Pattern.compile(PEAK_STR);

    final Peak peak;
    final Atom[] atoms;

    private DynamicsSource(Peak peak, Atom[] atoms) {
        this.atoms = atoms;
        this.peak = peak;
    }

    @Override
    public String toString() {
        StringBuilder sBuilder = new StringBuilder();
        sBuilder.append(peak.getName());
        for (Atom atom : atoms) {
            sBuilder.append(" ");
            sBuilder.append(atom.getFullName());
        }
        return sBuilder.toString();
    }

    public static Peak getPeak(String peakSpecifier, int nDim) {
        Matcher matcher = PEAK_PATTERN.matcher(peakSpecifier);
        if (!matcher.matches()) {
            throw new IllegalArgumentException("Invalid peak Specifier " + peakSpecifier);
        }
        PeakList peakList = PeakList.get(matcher.group(1));
        if (peakList == null) {
            peakList = new PeakList(peakSpecifier.substring(0, peakSpecifier.indexOf(".")), nDim);
        }
        int peakId = Integer.parseInt(matcher.group(2));
        Peak peak = peakList.getPeakByID(peakId);
        if (peak == null) {
            peak = peakList.getNewPeak();
            peak.setIdNum(peakId);
        }
        return peak;
    }

    public static Atom getAtom(String molName, String atomSpecifier, int peakID) {
        MoleculeBase molecule;

        if ((molName == null) || molName.equals("")) {
            molecule = MoleculeFactory.getActive();
        } else {
            molecule = MoleculeFactory.getMolecule(molName);
        }
        if (molecule == null) {
            molecule = MoleculeFactory.newMolecule(molName);
        }
        Matcher matcher = RESIDUE_PATTERN.matcher(atomSpecifier);
        String chainName = "";
        String resChar = "";
        String resNumStr = "";
        String ijChar = "";
        String atomName = "";
        if (matcher.matches()) {
            chainName = matcher.group(2);
            resChar = matcher.group(4);
            resNumStr = matcher.group(5);
            ijChar = matcher.group(6);
            atomName = matcher.group(8);
        }
        String resName = resChar == null ? "X" : PDBAtomParser.convert1To3(resChar);

        Polymer polymer = null;
        if (!molecule.getPolymers().isEmpty()) {
            if (chainName == null) {
                polymer = molecule.getPolymers().get(0);
            } else {
                polymer = molecule.getPolymer(chainName);
            }
        }
        if (polymer == null) {
            polymer = new Polymer(chainName);
            molecule.addEntity(polymer);
        }
        if (ijChar != null) {
            resNumStr = String.valueOf(peakID);
            resName = "X";
        }
        Residue residue = polymer.getResidue(resNumStr);
        if (residue == null) {
            residue = new Residue(resNumStr, resName);
            polymer.addResidue(residue);
        }
        Atom atom = residue.getAtom(atomName);
        if (atom == null) {
            atom = Atom.genAtomWithElement(atomName, atomName.substring(0, 1));
            residue.addAtom(atom);
        }
        return atom;
    }

    public static DynamicsSource createFromSpecifiers(String peakSpecifier,
            String resSpecifier, String... atomNames) {

        Peak peak = getPeak(peakSpecifier, atomNames.length);
        Atom[] atoms = new Atom[atomNames.length];
        String molName = null;
        int iAtom = 0;
        for (String atomName : atomNames) {
            Atom atom = getAtom(molName, resSpecifier + "." + atomName, iAtom);
            atoms[iAtom++] = atom;
        }
        return new DynamicsSource(peak, atoms);
    }

    public static DynamicsSource createFromPeak(Peak peak) {
        String molName = null;
        PeakList peakList = peak.getPeakList();
        int nDim = peakList.getNDim();
        Atom[] atoms = new Atom[nDim];
        for (int i = 0; i < nDim; i++) {
            String label = peak.getPeakDim(i).getAtomLabel();
            if ((label == null) || (label.trim().equals(""))) {
                label = peak.getPeakDim(i).getLabel();
                if ((label == null) || (label.trim().equals(""))) {
                    String fullPattern = peakList.getSpectralDim(i).getPattern();
                    if (fullPattern.length() > 2) {
                        String[] patterns = SpectralDim.parsePattern(fullPattern);
                        label = patterns[0];
                    } else {
                        String atomType = peak.getPeakDim(i).getSpectralDimObj().getAtomType();
                        label = "i." + atomType;
                    }
                }
            }
            Atom atom = getAtom(molName, label, peak.getIdNum());
            atoms[i] = atom;
        }
        return new DynamicsSource(peak, atoms);
    }

}

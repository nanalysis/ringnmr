/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.data;

import java.util.Optional;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.nmrfx.chemistry.Atom;
import org.nmrfx.chemistry.MoleculeBase;
import org.nmrfx.chemistry.MoleculeFactory;
import org.nmrfx.chemistry.Polymer;
import org.nmrfx.chemistry.Residue;
import org.nmrfx.chemistry.io.PDBAtomParser;
import org.nmrfx.chemistry.relax.ResonanceSource;
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

    final boolean createMol;
    final boolean createAtom;
    final boolean createPeakList;
    final boolean createPeak;

    public DynamicsSource(boolean createPeakList, boolean createPeak, boolean createMol, boolean createAtom) {
        this.createPeakList = createPeakList;
        this.createPeak = createPeak;
        this.createMol = createMol;
        this.createAtom = createAtom;
    }

    public Optional<Peak> getPeak(String peakSpecifier, int nDim) {
        Optional<Peak> result = Optional.empty();
        Matcher matcher = PEAK_PATTERN.matcher(peakSpecifier);
        if (!matcher.matches()) {
            throw new IllegalArgumentException("Invalid peak Specifier " + peakSpecifier);
        }
        PeakList peakList = PeakList.get(matcher.group(1));
        if (peakList == null) {
            if (createPeakList) {
                peakList = new PeakList(peakSpecifier.substring(0, peakSpecifier.indexOf(".")), nDim);
            } else {
                return result;
            }
        }
        int peakId = Integer.parseInt(matcher.group(2));
        Peak peak = peakList.getPeakByID(peakId);
        if (peak == null) {
            if (createPeak) {
                peak = peakList.getNewPeak();
                peak.setIdNum(peakId);
            } else {
                return result;
            }
        }
        return Optional.of(peak);
    }

    public Optional<Atom> getAtom(String molName, String atomSpecifier,
            int peakID, boolean patLabel) {
        MoleculeBase molecule;
        Optional<Atom> result = Optional.empty();

        if ((molName == null) || molName.equals("")) {
            molecule = MoleculeFactory.getActive();
        } else {
            molecule = MoleculeFactory.getMolecule(molName);
        }
        if (molecule == null) {
            if (!createMol) {
                return result;
            }
            if (molName == null) {
                molName = "noname";
            }
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
        return Optional.of(atom);
    }

    public Optional<ResonanceSource> createFromSpecifiers(String peakSpecifier,
            String resSpecifier, String... atomNames) {
        Optional<ResonanceSource> empty = Optional.empty();

        Optional<Peak> peakOpt = getPeak(peakSpecifier, atomNames.length);
        if (!peakOpt.isPresent()) {
            return empty;
        }
        Atom[] atoms = new Atom[atomNames.length];
        String molName = null;
        int iAtom = 0;
        for (String atomName : atomNames) {
            Optional<Atom> atomOpt = getAtom(molName, resSpecifier + "." + atomName, iAtom, false);
            if (atomOpt.isPresent()) {
                atoms[iAtom++] = atomOpt.get();
            } else {
                return empty;
            }
        }
        return Optional.of(new ResonanceSource(peakOpt.get(), atoms));
    }

    public Optional<ResonanceSource> createFromAtom(String peakSpecifier,
            Atom atom) {
        Optional<ResonanceSource> empty = Optional.empty();
        Atom[] atoms = {atom};
        Optional<Peak> peakOpt = getPeak(peakSpecifier, atoms.length);
        if (!peakOpt.isPresent()) {
            return empty;
        }
        return Optional.of(new ResonanceSource(peakOpt.get(), atoms));
    }

    public Optional<ResonanceSource> createFromAtomSpecifiers(String peakSpecifier,
            String... atomSpecifiers) {
        Optional<ResonanceSource> empty = Optional.empty();

        Optional<Peak> peakOpt = getPeak(peakSpecifier, atomSpecifiers.length);
        if (!peakOpt.isPresent()) {
            return empty;
        }
        Atom[] atoms = new Atom[atomSpecifiers.length];
        String molName = null;
        int iAtom = 0;
        for (String atomSpecifier : atomSpecifiers) {
            Optional<Atom> atomOpt = getAtom(molName, atomSpecifier, iAtom, false);
            if (atomOpt.isPresent()) {
                atoms[iAtom++] = atomOpt.get();
            } else {
                return empty;
            }
        }
        return Optional.of(new ResonanceSource(peakOpt.get(), atoms));
    }

    public Optional<ResonanceSource> createFromPeak(Peak peak) {
        Optional<ResonanceSource> empty = Optional.empty();
        String molName = null;
        PeakList peakList = peak.getPeakList();
        int nDim = peakList.getNDim();
        Atom[] atoms = new Atom[nDim];
        for (int i = 0; i < nDim; i++) {
            String label = peak.getPeakDim(i).getAtomLabel();
            boolean patLabel = false;
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
                    patLabel = true;
                }
            }
            Optional<Atom> atomOpt = getAtom(molName, label, peak.getIdNum(), patLabel);
            if (atomOpt.isPresent()) {
                atoms[i] = atomOpt.get();
            } else {
                return empty;
            }
        }
        return Optional.of(new ResonanceSource(peak, atoms));
    }

}

import numpy as np
from prody import *


class GammaResidue(Gamma):
    """Return force constant for a particular pair of residues.
    Returns normal force constant apart from at the residues being tested for allosteric character.
    """
    # see ProDy documentation for template
    def __init__(self, atoms, res_nos, frac_change, gamma=1):
        rnum = atoms.getResindices()
        # residue indices
        self._rnum = rnum
        self._res_nos = res_nos
        self._frac_change = float(frac_change)
        self._gamma = float(gamma)

    def gamma(self, dist2, i, j):
        rnum = self._rnum
        res_nos = self._res_nos
        gamma = self._gamma
        # if residue number is one of the input residues returns the modified force constant
        if rnum[i] in res_nos or rnum[j] in res_nos:
            return gamma * self._frac_change
        # otherwise return the normal force constant
        else:
            return gamma


def quantify_difference(anm_active, anm_mod_active, no_modes):
    """Compares standard to perturbed normal modes at the active site residues and returns a single perturbation value."""
    evals = anm_active.getEigvals()[:no_modes]
    # average magnitude difference for each normal mode
    mags_overall = []
    # iterate over low-frequency modes
    for mode_no in range(no_modes):
        # get mode and modified mode
        mode = anm_active[mode_no].getEigvec()
        mode_mod = anm_mod_active[mode_no].getEigvec()
        # sometimes the modified mode is reversed, so also compare the mode with the reversed modified mode
        mode_mod_alt = -anm_mod_active[mode_no].getEigvec()
        # add the minimum of the differences
        av_mags = vector_difference(mode, mode_mod)
        av_mags_alt = vector_difference(mode, mode_mod_alt)
        # add the minimum of the differences, which will correspond to the correct orientation
        site_av = min(av_mags, av_mags_alt)
        mags_overall.append(site_av)
    av = weight_by_eval(evals, mags_overall)
    return av


def vector_difference(mode1, mode2):
    """Returns average magnitude of vector difference between two eigenvectors."""
    mags = []
    # iterate over each atom
    for atom_index in range(0,len(mode1),3):
        # get normal and perturbed vectors
        vec1 = mode1[atom_index:atom_index+3]
        vec2 = mode2[atom_index:atom_index+3]
        # calculate magnitude of vector linking vectors
        diff = vec1 - vec2
        mag = np.sqrt(diff.dot(diff))
        mags.append(mag)
    return np.average(mags)


def weight_by_eval(evals, diffs):
    """Returns the average perturbation across the normal modes weighted by eigenvalue."""
    weightings = 1 / np.sqrt(evals)
    sum_weightings = sum(weightings)
    mags_weighted = weightings*diffs
    weighted_av = sum(mags_weighted) / sum_weightings
    return weighted_av


import numpy as np
from ..constants import *
from .design import calc_mz

# dummy value for presence/absence attributes
TAG = 'tag'

class DESIGN_6():
    name = 'pool1_termR'
    rule_set = 'pool0_termR'
    exclude_regex = 'NG'

    min_length = 9
    max_length = 13
    num_to_generate = int(1e6) # test
    num_generation_runs = 10 # test

    # does this stay the same?
    precursor_mz_min = 500.2667
    precursor_mz_max = 850
    precursor_bin_width = 10 # fewer mz bins for less prosit jobs
    precursor_bins = np.arange(precursor_mz_min, precursor_mz_max, precursor_bin_width)
    precursor_bin_names = {x: '{:.2f}'.format(x) for x in precursor_bins}

    # sometimes Prosit predicts iRT between 0 and 120 but the peptide elutes
    # between -25 and -10
    iRT_bin_width = 8
    iRT_bins = np.arange(-5, 150, iRT_bin_width)
    iRT_bin_names = {x: '{:.1f}'.format(x) for x in iRT_bins}
    normalized_collision_energy = 27

    # number of peptides retained for barcode selection, increase to saturate downstream filters
    input_barcodes_per_iRT_mz_bin = int(1e6)
    pred_barcodes_per_mz_bin = (
        input_barcodes_per_iRT_mz_bin * 
        len(iRT_bins) * 2)

    # filter peptides by number of product ions
    min_usable_ions = 5
    usable_ion_intensity = 0.25
    ion_mz_start = 200.1517
    ion_mz_max = 1300
    usable_ion_gate = (f'ion_type == "y" & {ion_mz_start} < ion_mz < {ion_mz_max}'
        f' & ion_charge == 1 & 2 < ion_length < (length - 1)'
        f' & {usable_ion_intensity} < intensity_prosit')
    # exact values don't matter for filtering by number of ions
    ignore_ion_intensity = 0.25
    ion_bins = np.arange(ion_mz_start, ion_mz_max, MZ_DOUBLE_SPACING)
    ion_bin_width = MZ_DOUBLE_SPACING

    ms1_resolution = 30000, 70000, 200000, int(1e6)


class DESIGN_7(DESIGN_6):
    name = 'pool1_termK'
    rule_set = 'pool0_termK'

class DESIGN_8(DESIGN_6):
    name = 'chip162_termR'
    rule_set = 'chip162_termR'
    min_length = 6
    max_length = 13


class DESIGN_9(DESIGN_6):
    name = 'chip162_termK'
    rule_set = 'chip162_termK'
    min_length = 6
    max_length = 13


class DESIGN_10(DESIGN_6):
    """Exhaustively sample short barcodes
    """
    name = 'chip162_termK'
    rule_set = 'chip162_termK'
    exhaustive = TAG
    min_length = 9
    max_length = 9
    precursor_mz_min = 450
    precursor_mz_max = 600
    precursor_bin_width = 6 # fewer mz bins for less prosit jobs
    pred_barcodes_per_mz_bin = int(1e5)
    peptide_gate = 'hydro_levy2012 < 2.5'


class DESIGN_11(DESIGN_6):
    """Standard Nterm barcodes with new rule (no Y, L=>I) and up-front hydrophobicity filtering.
    """
    name = 'chip162_termK'
    rule_set = 'chip162_termK'
    min_length = 9
    max_length = 13
    peptide_gate = 'hydro_levy2012 < 2.5'

class DESIGN_12(DESIGN_10):
    name = 'chip162_termR_short'
    rule_set = 'chip162_termR'

class DESIGN_13(DESIGN_11):
    """C-term barcodes similar to DESIGN_11.
    """
    name = 'chip162_termR'
    rule_set = 'chip162_termR'
    peptide_gate = 'hydro_levy2012 < 2.5'



class DESIGN_14(DESIGN_6):
    """Exhaustively sample short Nterm barcodes
    """
    name = 'chip162_termK'
    rule_set = 'chip162_termK'
    exhaustive = TAG
    min_length = 8
    max_length = 8
    precursor_mz_min = 450
    precursor_mz_max = 600
    precursor_bin_width = 6 # fewer mz bins for less prosit jobs
    pred_barcodes_per_mz_bin = int(1e5)
    peptide_gate = 'hydro_levy2012 < 2.5'


class DESIGN_TEST(DESIGN_6):
    num_to_generate = 1000
    num_generation_runs = 1
    # does this stay the same?
    precursor_mz_min = 500.2667
    precursor_mz_max = 550
    precursor_bin_width = 10  # fewer mz bins for less prosit jobs
    precursor_bins = np.arange(
        precursor_mz_min, precursor_mz_max, precursor_bin_width)
    precursor_bin_names = {x: '{:.2f}'.format(x) for x in precursor_bins}

import numpy as np
from .constants import *
from . import flycodes


class DESIGN_0():
    name = 'random-termK'
    min_length = 10
    max_length = 13
    num_to_generate = int(4e6)
    num_generation_runs = 30

    precursors_per_bin = 200
    precursor_mz_start = 500.2667
    precursor_mz_max = 850
    precursor_bins = np.arange(precursor_mz_start, precursor_mz_max, MZ_DOUBLE_SPACING)
    precursor_bin_width = MZ_DOUBLE_SPACING
    precursor_bin_names = {x: '{:.2f}'.format(x) for x in precursor_bins}

    iRT_bins = np.arange(-10, 110, 5)
    iRT_bin_names = {x: '{:.1f}'.format(x) for x in iRT_bins}
    iRT_bin_width = 5

    normalized_collision_energy = 27

    # number of peptides retained for barcode selection
    # tested empirically to saturate barcode selection
    input_barcodes_per_iRT_mz_bin = 10000

    ion_mz_start = 200.1517
    ion_mz_max = 1300
    usable_ion_intensity = 0.1
    ignore_ion_intensity = 0.01
    usable_ion_gate = ('ion_type == "y" & {} < ion_mz < {} & ion_charge == 1'
        '& 2 < ion_length < (length - 1) & {} < intensity_prosit'
        .format(ion_mz_start, ion_mz_max, usable_ion_intensity))

    ion_bins = np.arange(ion_mz_start, ion_mz_max, MZ_DOUBLE_SPACING)
    ion_bin_width = MZ_DOUBLE_SPACING

    selection_seeds = np.arange(15)
    min_unique_ions = 2



class DESIGN_1():
    name = 'RJ_76'

    precursor_bins = np.linspace(550, 850, 100)
    precursor_bin_width = np.diff(precursor_bins)[0] # no space between bins
    precursor_bin_min_spacing = 0
    precursors_per_bin = 200

    iRT_bins = np.linspace(-25, 190, 11) # wider than the random peptide iRT range
    iRT_bin_width = np.diff(iRT_bins)[0] # no space between bins
    iRT_bin_min_spacing = 0

    if np.diff(precursor_bins).min() < precursor_bin_min_spacing:
        raise ValueError









import numpy as np

class DESIGN_0():
    name = 'random_peptides'
    min_length = 10
    max_length = 13
    num_to_generate = int(1e6)
    num_permutations = 300
    # min_spacing = 0.15

    precursor_bins = np.linspace(550, 850, 100)
    precursor_bin_width = 1
    precursor_bin_min_spacing = 2.5
    precursors_per_bin = 200

    iRT_bins = np.linspace(-10, 110, 11)
    iRT_bin_width = 6
    iRT_bin_min_spacing = 6

    if np.diff(precursor_bins).min() < precursor_bin_min_spacing:
        raise ValueError


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



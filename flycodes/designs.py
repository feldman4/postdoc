import numpy as np
from ..constants import *
from .design import calc_mz


class DESIGN_0():
    """First large design, produces barcodes for every iRT,mz bin that can be
    downsampled (e.g., to avoid isotope overlap) later.
    """
    name = 'random-termK'
    rule_set = 'RJ_noH_termK'
    min_length = 10
    max_length = 13
    num_to_generate = int(4e6)
    num_generation_runs = 30

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
    input_barcodes_per_iRT_mz_bin = 8000
    pred_barcodes_per_mz_bin = (
        input_barcodes_per_iRT_mz_bin * 
        len(iRT_bins) * 2)
    
    ion_mz_start = 200.1517
    ion_mz_max = 1300
    usable_ion_intensity = 0.1
    ignore_ion_intensity = 0.01
    usable_ion_gate = ('ion_type == "y" & {} < ion_mz < {} & ion_charge == 1'
        '& 2 < ion_length < (length - 1) & {} < intensity_prosit'
        .format(ion_mz_start, ion_mz_max, usable_ion_intensity))

    ion_bins = np.arange(ion_mz_start, ion_mz_max, MZ_DOUBLE_SPACING)
    ion_bin_width = MZ_DOUBLE_SPACING

    selection_seeds = np.arange(5)
    min_unique_ions = 2

    exclude_regex = 'NG'


class DESIGN_1(DESIGN_0):
    """C-term barcodes.
    """
    name = 'random-termR'
    rule_set = 'RJ_noH_termR'


class DESIGN_2(DESIGN_0):
    """Based on permutations of known good peptides.
    """
    name = 'RJ_76_termR'
    rule_set = 'RJ_76'

    num_generation_runs = 1
    num_permutations = 1000
    known_peptides = ['AGNLGGGVVTLER', 'ALLAAQYSGAQVR', 'ALQNAVTTFVNR', 
    'APLLLATDVASR', 'AQFEGLVTDLLR', 'AQLFALTGVQPAR', 'AQLFANTVDNAR', 
    'AQLLQPTLELNPR', 'ASGQAFELLLSPR', 'AVLGEVVLYSGAR', 'DAGQLSGLNVLR', 
    'DAGTLAGLNVLR', 'DAGVLAGLNVLR', 'DAPEEEDHVLVLR', 'DAVVYPLLVEFTR', 
    'DFSPSGLFGAFQR', 'DGETPDPEDPSR', 'DSALEFLTQLSR', 'DSAQTSVTQAQR', 
    'DTQSGSLLFLGR', 'EDLTQSAQHALR', 'EVDLGLPDATGR', 'EVTLNQSLLAPLR', 
    'EVYQQQQYGSGGR', 'FAAATGATPLAGR', 'FDDGAGGDNEVQR', 'FEELNADLFR', 
    'FGSDQSENVDR', 'GATQQLLDEAER', 'GLTPSQLGVLLR', 'GPQVQQPPPSNR', 
    'GQLEALQVDGGR', 'GVQVETLSPGDGR', 'LASGLGLAWLVGR', 'LATLLGLQAPPTR', 
    'LDLDSPPLTAR', 'LDYLGVSYGLTPR', 'LEGDETSTEAATR', 'LFVTNDAATLLR', 
    'LGEYGFQNALLVR', 'LLDTLGLSQPQWR', 'LLTSFLPAQLLR', 'LPGLLLAASAVR', 
    'LPVGTTATLYFR', 'LPWFQYPLLYDLR', 'LQSSQEPEAPPPR', 'LQTQPGYANTLR', 
    'LSGLLYEETR', 'LSVNSVTAGDYSR', 'LVALVDVLDQNR', 'LVDQNLFSFYLSR', 
    'NLFPSNLVSAAFR', 'NPAVLSAASFDGR', 'NPQNSSQSADGLR', 'NVNSNLLTWNR', 
    'SAYGGPVGAGLR', 'SELGNQSPSTSSR', 'SFPDFPTPGVVFR', 'SLGSVQAPSYGAR', 
    'SNVSDAVAQSTR', 'SQLLQYVYNLVPR', 'TLSFGSDLNYATR', 'TNEAQALETAR', 
    'TPFLLVGTQLDLR', 'TTPDVLFVFGFR', 'TYLNPFVSFLDQR', 'VALTGLTVAEYFR', 
    'VEPGLGADNSVVR', 'VLANPGNSQVAR', 'VQLSPDSGGLPER', 'WSNLPFLTVPLSR',
    'WVAVVVPSGQEQR','YGVNPGPLVGTTR', 'YNSQLLSFVR', 'YTLSQEAYDQR', 
    'YYVTLLDAPGHR']

    PARENT = DESIGN_0
    mz = [calc_mz(x, 2) for x in known_peptides]
    mask = (np.abs(PARENT.precursor_bins[:, None] - mz) 
                   < PARENT.precursor_bin_width).any(axis=1)
    precursor_bins = PARENT.precursor_bins[mask]
    precursor_bin_names = {x: '{:.2f}'.format(x) for x in precursor_bins}


class DESIGN_3(DESIGN_0):
    parent = DESIGN_0
    name = 'pool0_termK'
    rule_set = 'pool0_termK'

    iRT_bins = np.arange(-25, 150, 5)
    iRT_bin_names = {x: '{:.1f}'.format(x) for x in iRT_bins}
    iRT_bin_width = 5

    # MS1 selection doesn't benefit much from re-runs
    selection_seeds = [0]

    # subdivide MS1 bins into N groups
    ms1_selection_scans = 10

    n = int(np.ceil(len(parent.precursor_bin_names) / ms1_selection_scans))
    ms1_selection_ranges = {}
    for i in range(ms1_selection_scans):
        bins = list(parent.precursor_bin_names.values())[i * n:(i + 1) * n]
        key = bins[0] + '-' + bins[-1]
        ms1_selection_ranges[key] = bins
    del n, i, bins, key

    # set high enough to ensure barcodes found for all precursor masses
    ms1_selection_input_max = 200


class DESIGN_4(DESIGN_3):
    name = 'pool0_termR'
    rule_set = 'pool0_termR'

    # num_to_generate = 1e4
    # num_generation_runs = 10

    # precursor_bins = np.arange(
    #     DESIGN_0.precursor_mz_start, 
    #     DESIGN_0.precursor_mz_max, 
    #     MZ_DOUBLE_SPACING
    #     )[100:101]
    # precursor_bin_names = {x: '{:.2f}'.format(x) for x in precursor_bins}

    # iRT_bins = np.arange(-10, 110, 5)[10:11]
    # iRT_bin_names = {x: '{:.1f}'.format(x) for x in iRT_bins}


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
    usable_ion_gate = ('ion_type == "y" & {} < ion_mz < {} & ion_charge == 1'
        '& 2 < ion_length < (length - 1) & {} < intensity_prosit'
        .format(ion_mz_start, ion_mz_max, usable_ion_intensity))
    # exact values don't matter for filtering by number of ions
    ignore_ion_intensity = 0.25
    ion_bins = np.arange(ion_mz_start, ion_mz_max, MZ_DOUBLE_SPACING)
    ion_bin_width = MZ_DOUBLE_SPACING

    ms1_resolution = 30000, 70000, 200000, int(1e6)


class DESIGN_7(DESIGN_6):
    name = 'pool1_termK'
    rule_set = 'pool0_termK'


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


class DESIGN_WW():
    name = 'pool1_terK'
    rule_set = 'pool0_termK'
    exclude_regex = 'NG'

    min_length = 6
    max_length = 6
    num_to_generate = int(1e3) # test
    num_generation_runs = 10 # test

    # does this stay the same?
    precursor_mz_min = 300#500.2667
    precursor_mz_max = 400
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
    usable_ion_gate = ('ion_type == "y" & {} < ion_mz < {} & ion_charge == 1'
        '& 2 < ion_length < (length - 1) & {} < intensity_prosit'
        .format(ion_mz_start, ion_mz_max, usable_ion_intensity))
    # exact values don't matter for filtering by number of ions
    ignore_ion_intensity = 0.25
    ion_bins = np.arange(ion_mz_start, ion_mz_max, MZ_DOUBLE_SPACING)
    ion_bin_width = MZ_DOUBLE_SPACING

    ms1_resolution = 30000, 70000, 200000, int(1e6)


runs = {
    'run_001': DESIGN_0,
    'run_002': DESIGN_1,
    'run_003': DESIGN_3,
    'run_004': DESIGN_4,
    'run_005': DESIGN_2,
    'run_006': DESIGN_6,
    'run_007': DESIGN_6,
    'run_008': DESIGN_7,
    'run_test': DESIGN_TEST,
    'run_ww': DESIGN_WW,
    }


# fail dnachisel constraints
exclude_from_synthesis = [
    'PWGVQQVVVSEAR',
    'QSVVQAWGEVPVR',
    'SAEQQWGVVVVPR',
    'QVVPEVSAWGQVR',
    'ESAQVQPWGVVVR',
    'DYFTSSTAGLLNR',
    'STGSAFNLDYTLR',
    'DYSFTGASLTLNR',
    'STYGSFTDLNLAR',
    'SFTLNSLDGAYTR',
    ]

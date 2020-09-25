import pyteomics.ms1
import pyteomics.mzml
import pandas as pd
import numpy as np
import os
import gzip
import re

from ..constants import HOME

class Ecoli:
    def __init__(self):
        with gzip.open('flycodes/ecoli_ref/BL21.gff3', 'rb') as fh:
            self.gff3 = fh.read().decode()

    def download_ref():
        base = ('ftp://ftp.ensemblgenomes.org/pub/bacteria/release-47/'
                '{suffix}/bacteria_11_collection/escherichia_coli_bl21_de3_/')
        base_gff3 = base.format(suffix='gff3')
        base_fasta = base.format(suffix='fasta')
        release = {
            'gff3': f'{base_gff3}Escherichia_coli_bl21_de3_.ASM956v1.46.gff3.gz',
            'pep.fa': f'{base_fasta}pep/Escherichia_coli_bl21_de3_.ASM956v1.pep.all.fa.gz'
        }

        for suffix, url in release.items():
            local = f'flycodes/ecoli_ref/BL21.{suffix}.gz'
            try:
                os.remove(local)
            except OSError:
                pass
            wget.download(url, local)

    def find_gene_symbol(self, symbol):
        pat = f'id=gene:(\w+);name={symbol.lower()}'
        return re.findall(pat, self.gff3.lower())[0].upper()

    def get_peptide_sequences(self, symbols):
        from ..sequence import read_fasta

        records = read_fasta('flycodes/ecoli_ref/BL21.pep.fa.gz')
        gene_ids = {}
        for s in symbols:
            gene_ids[s] = self.find_gene_symbol(s)

        sequences = []
        for k, v in gene_ids.items():
            for name, seq in records:
                if v in name:
                    sequences += [(k, seq)]

        return (pd.DataFrame(sequences)
                .rename(columns={0: 'gene_symbol', 1: 'seq_aa'})
                )

    def get_crap_table(self, drive):
        return (drive('mass spec barcoding/contaminants')
            .pipe(lambda x:
                    x.merge(self.get_peptide_sequences(
                        x['gene_symbol'].dropna())))
            [['gene_symbol', 'seq_aa']]
            .drop_duplicates()
            )


def e_coli_crap_mz(drive, metadata):
    df_crap = Ecoli().get_crap_table(drive)

    peptides = set()
    for s in df_crap['seq_aa']:
        s_ = s.replace('R', 'R.').replace('K', 'K.')
        peptides |= set(s_.split('.'))
    peptides -= {''}

    mz_list = []
    for p in peptides:
        mz_list.append(calc_mz(p, charge=2))

    bins = bin_by_value(mz_list, metadata.precursor_bins, metadata.precursor_bin_width)
    bad_bins = pd.Series(bins).dropna().drop_duplicates().pipe(sorted)
    return bad_bins



def load_mzml_data(filename, progress=lambda x: x):
    mz = []
    intensity = []
    info = []
    reader = pyteomics.mzml.MzML(filename)
    for spectrum in progress(reader):
        scan = spectrum['scanList']['scan'][0]
        start = scan['scan start time']
        inject = scan['ion injection time']
        peak_
        mz_arr = spectrum['m/z array']
        int_arr = spectrum['intensity array']
        info += [{'scan': start, 'inject': inject}]
        mz += [mz_arr]
        intensity += [int_arr]

    mz = np.array(mz)
    intensity = np.array(intensity)
    df_info = pd.DataFrame(info)
    return mz, intensity, df_info

def load_ms1_data(filename, progress=lambda x: x):
    mz = []
    intensity = []
    info = []
    reader = pyteomics.ms1.MS1(filename)
    keys = ['RTime', 'BPI', 'BPM', 'TIC']
    for scan in progress(reader):
        mz_arr = scan['m/z array']
        int_arr = scan['intensity array']
        info += [scan['params']]
        mz += [mz_arr]
        intensity += [int_arr]

    mz = np.array(mz)
    intensity = np.array(intensity)
    df_info = pd.DataFrame(info)
    return mz, intensity, df_info

def grid_ms1_intensities(mz, intensity, time):
    """Rounds down mz and time. Multiply input and divide output to increase precision.
    Integrate intensity falling into each time,mz bin.
    """
    time = np.array(time).astype(int)

    all_mz = np.concatenate([x.flatten() for x in mz])
    all_int = np.concatenate([x.flatten() for x in intensity])

    mz_arr = np.zeros((time[-1] + 1, int(all_mz.max()) + 1))
    for mz_, int_, rt in zip(mz, intensity, time):
        mz_arr[rt, mz_.astype(int)] += int_
    
    mz_arr[mz_arr == 0] = np.nan
    
    df_int = pd.DataFrame(mz_arr).dropna(axis=1, how='all')
    df_int.index.name = 'time'
    df_int.columns.name = 'mz'
    return df_int


def generate_msfragger_cmd(mzML, protein_fa, output_format='pepXML'):
    """
    output_format can be tsv or pepXML
    """
    params_text = msfragger_params.format(
        protein_fa=protein_fa, output_format=output_format)
    params = HOME / 'Downloads' / 'fragger.test.params'
    with open(params, 'w') as fh:
        fh.write(params_text)
    
    java_flags = '-jar -Dfile.encoding=UTF-8 -Xmx6G'
    cmd = f'java {java_flags} {msfragger_jar} {params} {mzML}'
    return cmd


msfragger_jar = 'misc/msfragger/MSFragger-3.0.jar'
msfragger_params = """database_name = {protein_fa}
num_threads = 3			# 0=poll CPU to set num threads; else specify num threads directly (max 64)

precursor_mass_lower = -10
precursor_mass_upper = 10
precursor_mass_units = 1			# 0=Daltons, 1=ppm
precursor_true_tolerance = 10
precursor_true_units = 1			# 0=Daltons, 1=ppm
fragment_mass_tolerance = 10
fragment_mass_units = 1			# 0=Daltons, 1=ppm
calibrate_mass = 0			# 0=Off, 1=On, 2=On and find optimal parameters
decoy_prefix = rev_

deisotope = 1
isotope_error = -1/0/1/2			# 0=off, -1/0/1/2/3 (standard C13 error)
mass_offsets = 0			# allow for additional precursor mass window shifts. Multiplexed with isotope_error. mass_offsets = 0/79.966 can be used as a restricted ‘open’ search that looks for unmodified and phosphorylated peptides (on any residue)
precursor_mass_mode = selected

remove_precursor_peak = 0			# 0 = not remove, 1 = only remove the peak with the precursor charge, 2 = remove all peaks with all charge states. Default: 0
remove_precursor_range = -1.50,1.50			# Unit: Da. Default: -1.5,1.5
intensity_transform = 0			# 0 = none, 1 = sqrt root. Default: 0

localize_delta_mass = 0
delta_mass_exclude_ranges = (-1.5,3.5)
fragment_ion_series = b,y

search_enzyme_name = trypsin
search_enzyme_cutafter = KR
search_enzyme_butnotafter = P

num_enzyme_termini = 2			# 2 for enzymatic, 1 for semi-enzymatic, 0 for nonspecific digestion
allowed_missed_cleavage = 1			# maximum value is 5

clip_nTerm_M = 0

#maximum of 7 mods - amino acid codes, * for any amino acid, [ and ] specifies protein termini, n and c specifies peptide termini
variable_mod_01 = 15.99490 M 3
variable_mod_02 = 42.01060 [^ 1
# variable_mod_03 = 79.96633 STY 3
# variable_mod_04 = -17.02650 nQnC 1
# variable_mod_05 = -18.01060 nE 1
# variable_mod_06 = 0.00000 site_06 1
# variable_mod_07 = 0.00000 site_07 1

allow_multiple_variable_mods_on_residue = 0			# static mods are not considered
max_variable_mods_per_mod = 3
max_variable_mods_per_peptide = 3			# maximum 5
max_variable_mods_combinations = 5000			# maximum 65534, limits number of modified peptides generated from sequence

output_file_extension = {output_format}         # pepXML or tsv
output_format = {output_format}
output_report_topN = 1
output_max_expect = 50
report_alternative_proteins = 0			# 0=no, 1=yes

precursor_charge = 1 4			# precursor charge range to analyze; does not override any existing charge; 0 as 1st entry ignores parameter
override_charge = 0			# 0=no, 1=yes to override existing precursor charge states with precursor_charge parameter

digest_min_length = 7
digest_max_length = 50
digest_mass_range = 200.0 5000.0			# MH+ peptide mass range to analyze
max_fragment_charge = 2			# set maximum fragment charge state to analyze (allowed max 5)

#open search parameters
track_zero_topN = 0			# in addition to topN results, keep track of top results in zero bin
zero_bin_accept_expect = 0.00			# boost top zero bin entry to top if it has expect under 0.01 - set to 0 to disable
zero_bin_mult_expect = 1.00			# disabled if above passes - multiply expect of zero bin for ordering purposes (does not affect reported expect)
add_topN_complementary = 0

# spectral processing

minimum_peaks = 10			# required minimum number of peaks in spectrum to search (default 10)
use_topN_peaks = 100
min_fragments_modelling = 2
min_matched_fragments = 4
minimum_ratio = 0.01			# filter peaks below this fraction of strongest peak
clear_mz_range = 0.0 0.0			# for iTRAQ/TMT type data; will clear out all peaks in the specified m/z range

# additional modifications

add_Cterm_peptide = 0.000000
add_Nterm_peptide = 0.000000
add_Cterm_protein = 0.000000
add_Nterm_protein = 0.000000

add_G_glycine = 0.000000
add_A_alanine = 0.000000
add_S_serine = 0.000000
add_P_proline = 0.000000
add_V_valine = 0.000000
add_T_threonine = 0.000000
add_C_cysteine = 57.021464
add_L_leucine = 0.000000
add_I_isoleucine = 0.000000
add_N_asparagine = 0.000000
add_D_aspartic_acid = 0.000000
add_Q_glutamine = 0.000000
add_K_lysine = 0.000000
add_E_glutamic_acid = 0.000000
add_M_methionine = 0.000000
add_H_histidine = 0.000000
add_F_phenylalanine = 0.000000
add_R_arginine = 0.000000
add_Y_tyrosine = 0.000000
add_W_tryptophan = 0.000000
add_B_user_amino_acid = 0.000000
add_J_user_amino_acid = 0.000000
add_O_user_amino_acid = 0.000000
add_U_user_amino_acid = 0.000000
add_X_user_amino_acid = 0.000000
add_Z_user_amino_acid = 0.000000
"""


import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from ..sequence import reverse_complement
from . import designs


def shifted_diagonals(df_ms):
    num_ms1, num_iRT = df_ms.shape
    df_ms = df_ms.copy()
    i = np.arange(num_ms1)
    
    # first diagonal
    start = np.arange(9) * int(num_iRT*1.1)
    for col in range(0, num_iRT, 2):
        df_ms.values[start + col, col] = f'scan_{col}'
    
    # second diagonal
    start = np.arange(9) * int(num_iRT*1.1)
    for col in range(1, num_iRT, 2):
        df_ms.values[start + col + int(num_iRT/2), col] = f'scan_{col}'
    
    return df_ms


def unshifted_diagonals(iRT, mz):
    df_ms = pd.DataFrame(index=mz, columns=iRT).fillna(False)

    num_ms1, num_iRT = df_ms.shape
    df_ms = df_ms.copy()
    iRT = 0
    for ms1 in range(num_ms1):
        iRT = (iRT + 2) % num_iRT
        df_ms.values[ms1, iRT] = True
    return df_ms.T


def limit_and_order_peptides(df_peptides, num_barcodes, bin_cols):
    """Request `num_barcodes` rows from `df_peptides`, keeping only as
    many bins as needed to get to `num_barcodes`. The order of rows
    in the returned table is obtained by cycling through sorted bins 
    one row at a time.
    """
    bins = df_peptides.groupby(bin_cols).size()
    bins = bins.sample(frac=1, random_state=0).rename('bin_size')
    up_to = np.where(bins.cumsum() > num_barcodes)[0][0]
    bins = bins.iloc[:up_to + 1]

    barcodes = (df_peptides
     .set_index(bin_cols)
     .join(bins, how='right')
    )

    barcodes['rank'] = (barcodes.assign(dummy=1)
     .groupby(bin_cols)['dummy'].rank(method='first'))

    barcodes = barcodes.sort_values(['rank'] + bin_cols)
    return barcodes.reset_index().drop('rank', axis=1)


def reverse_translate_barcode(seq_aa, gc_min=0.3, gc_max=0.7):
    sequence = dnachisel.biotools.reverse_translate(seq_aa)
    problem = dnachisel.DnaOptimizationProblem(
        sequence=sequence,
        constraints=[
            dnachisel.AvoidPattern("BsaI_site"),
            dnachisel.AvoidPattern('GGGG'),
            dnachisel.EnforceGCContent(mini=gc_min, maxi=gc_max),
            dnachisel.EnforceTranslation(),
        ],
        objectives=[dnachisel.CodonOptimize(species='e_coli')]
    )

    problem.resolve_constraints()
    assert problem.all_constraints_pass()
    return problem.sequence


def reverse_translate_barcodes(peptides, progress=lambda x: x):
    f = 'flycodes/pool0/dnachisel_output.csv'

    aa_to_dna = (pd.read_csv(f)
     .set_index('sequence_aa')['sequence_dna']
     .to_dict()            
    )

    bad_aa = []
    for x in progress(peptides):
        if x in aa_to_dna:
            continue
        try:
            aa_to_dna[x] = reverse_translate_barcode(x) 
        except:
            bad_aa += [x]
            
    return aa_to_dna


def build_oligos(df_oligos, oligo_templates, dialout_adapters):
    it = df_oligos[['vector', 'sequence_dna', 'subpool']].values
    arr = []
    for vector, seq, subpool in it:
        parts = oligo_templates[vector].split('N')
        cloning_fwd, cloning_rev = parts[0], parts[-1]
        pcr_fwd, pcr_rev = dialout_adapters[subpool]
        pcr_rev_rc = reverse_complement(pcr_rev)
        arr += [pcr_fwd + cloning_fwd + seq + cloning_rev + pcr_rev_rc]

    return arr


def select_pool0_peptides(df_layout):

    bad_peptides = designs.exclude_from_synthesis
    it = df_layout.groupby(['design', 'run', 'selection', 'mask'])

    arr = []
    for (design, run, selection, mask), df_layout_ in it:
        print(design, run, selection, mask, sep=', ')
    #     continue

        ms1 = '_ms1' in run
        run = run.replace('_ms1', '') if ms1 else run
        f = f'flycodes/{run}/barcode_ions.csv'        
        f = f.replace('.csv', '_ms1.csv') if ms1 else f
        mz_key = 'ms1_range' if ms1 else 'mz_bin'

        assert design == designs.runs[run].name


        # load table of ions used to identify barcodes
        # filter to one row per peptide
        num_barcodes = df_layout_['num_barcodes'].sum()
        df_barcodes = (pd.read_csv(f)
         .query('sequence != @bad_peptides')
         .drop_duplicates('sequence')
         .filter(regex='^(?!ion).*')
        )

        # apply mask to iRT,MS1 bins if provided
        if mask.endswith('csv'):
            df_ms = pd.read_csv(mask, index_col=0)
            df_ms.index.name = 'iRT_bin'
            df_ms.columns = (pd.Index(df_ms.columns, name='mz_bin')
                             .astype(float))
            masked_bins = df_ms.stack().loc[lambda x: x].rename('keep')

            df_barcodes = df_barcodes.join(
                masked_bins, on=['iRT_bin', 'mz_bin'], how='inner')

        # sort rows so that the first N/n barcodes belong in subpool 1 of n
        if selection == 'per iRT,MS1':
            df_barcodes = df_barcodes.pipe(
                limit_and_order_peptides, 
                num_barcodes, ['iRT_bin', mz_key])

        total_required = df_layout_['num_barcodes'].sum()
        assert total_required <= len(df_barcodes)
        df_barcodes = (df_barcodes[:total_required]
         .assign(subpool=np.repeat(df_layout_['subpool'], 
                                   df_layout_['num_barcodes']).values)
         .drop(['keep', 'bin_size', 'intensity_prosit'], axis=1, errors='ignore')
         .rename(columns={'run': 'snakemake_run'})
        )

        arr += [df_layout_.merge(df_barcodes)]
        
    df_oligos = (pd.concat(arr)
     .sort_values(['subpool', 'iRT_bin', 'mz_bin'])
    )
    
    return df_oligos
    

def build_ms():
    DESIGN = designs.DESIGN_3
    df_ms = (unshifted_diagonals(
        DESIGN.iRT_bins[1::2],
        DESIGN.precursor_bins)
    .rename(columns=DESIGN.precursor_bin_names)
    .reindex(index=DESIGN.iRT_bins).fillna(False)
    )
    df_ms.to_csv('flycodes/pool0/diagonal_iRT_MS1.csv')


def design_pool(drive):
    dialout_adapters = (drive('reagents/dialout_oligos')
     .dropna(subset=['dialout'])
     [['dialout', 'FWD_priming', 'REV_priming']]
     .pipe(lambda x: {int(y[0]): (y[1], y[2]) for y in x.values})
    )

    df_layout = drive('mass spec barcoding/pool0', skiprows=1)

    oligo_templates = (drive('reagents/plasmids')
     .set_index('name')['design_template'].to_dict())

    df_peptides = select_pool0_peptides(df_layout)

    aa_to_dna = reverse_translate_barcodes(df_peptides['sequence'])
    df_oligos = (df_peptides
    .assign(sequence_dna=lambda x: x['sequence'].map(aa_to_dna))
    .dropna(subset=['sequence_dna'])            
    .assign(oligo=lambda x: 
            x.pipe(build_oligos, oligo_templates, 
                    dialout_adapters))
    )

    return df_oligos


def plot_subpools_iRT_mz(df_oligos, progress=lambda x: x):
    height, width = 5, 9

    iRT_lim = df_oligos['iRT_bin'].describe()[['min', 'max']]
    mz_lim = df_oligos['mz_bin'].describe()[['min', 'max']]

    def plot(df):
        with sns.axes_style('whitegrid'):
            return (df.groupby(['iRT_bin', 'mz_bin']).size()
             .rename('count').reset_index()
             .pipe(sns.FacetGrid, hue='count', palette='Greens', height=height,
                  aspect=width/height)
             .map(plt.scatter, 'mz_bin', 'iRT_bin')
             .add_legend()
            )

    for _, df in progress(df_oligos.groupby('subpool')):
        fg = plot(df)
        ranges = df['ms1_range'].dropna().drop_duplicates()
        for i, ms1_range in enumerate(sorted(ranges)):
            left, right = ms1_range.split('-')
            left, right = float(left), float(right)
            color = [(0.9, 0.9, 0.9), (0.8, 0.8, 0.8)][i % 2]
            fg.ax.fill([left, right, right, left],
                [iRT_lim[0], iRT_lim[0], iRT_lim[1], iRT_lim[1]],
                color=color, fill=True, zorder=-1,
               )

        fg.ax.set_xlim(mz_lim)
        fg.ax.set_ylim(iRT_lim)
        title = 'subpool={subpool}; design={design}\n{description}'
        title = title.format(**df.iloc[0])
        title += f'\n{df.shape[0]}/{df.iloc[0]["num_barcodes"]} barcodes'
        fg.ax.set_title(title)

        f = 'flycodes/pool0/figures/pool0_{subpool}_iRT_mz_summary.png'
        f = f.format(**df.iloc[0])
        fg.savefig(f)
        plt.close(fg.fig)
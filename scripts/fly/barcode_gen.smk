"""
# example digs submission
sbatch -p gpu --mem=80g --gres=gpu:rtx2080:1 -c 10 s/fly/run.sh run_003
sbatch -p medium --mem=80g -c 10 s/fly/run.sh run_003

rm fly_runs.list
echo ~/s/fly/run.sh run_009 >> fly_runs.list
echo ~/s/fly/run.sh run_010 >> fly_runs.list
# app.sh submit fly_runs.list --queue=gpu --memory=80g -cpus=10 
"""

# IMPORTS
import postdoc.flycodes as fly
from postdoc.flycodes import designs
from postdoc.utils import timestamp, csv_frame

import numpy as np
import pandas as pd
import inspect

# CONSTANTS

GPU_MEM_FRACTION = 0.1

get_first = lambda xs: [x for x in xs if os.path.exists(x)][0]
MODEL_IRT = get_first([
    '/root/model_spectra', # apptainer
    '/projects/ms/prosit/models/irt/',
])
MODEL_SPECTRA = get_first([
    '/root/model_irt', # apptainer
    '/projects/ms/prosit/models/spectra/',
])

# CONFIG

M = designs.runs[config['run']]
RUN_NAME = f'{M.name}_{config["timestamp"]}' # unique name for this snakemake run
if 'exhaustive' in dir(M):
    RUNS = fly.design.get_prefixes(M.rule_set)
else:
    RUNS = ['{:03d}'.format(x) for x in range(M.num_generation_runs)]

def expand_ms1_range(wildcards):
    return expand('process/{{design}}_iRT_{{bin_iRT}}_mz_{bin_mz}.precursors.csv', 
        bin_mz=M.ms1_selection_ranges[wildcards.ms1_range])

rule all:
    input: 
        # get all barcodes, no MS1 filtering
        # 'all_barcodes.csv',

        # get barcodes at each MS1 resolution specified in design
        expand('barcodes_ms1_{ms1_res}.csv', ms1_res=M.ms1_resolution),

        # expand('process/{design}_{run}_{bin_mz}.peptides.csv', 
        #     design=M.name,
        #     run=RUNS,
        #     bin_mz=M.precursor_bin_names.values()),
        # expand('process/{design}_iRT_{bin_iRT}_mz_{bin_mz}.precursors.csv',
        #     design=M.name,
        #     bin_iRT=M.iRT_bin_names.values(),
        #     bin_mz=M.precursor_bin_names.values()),
        # expand('process/{design}_iRT_{bin_iRT}.ms1_{ms1_res}.csv',
        #     design=M.name,
        #     bin_iRT=M.iRT_bin_names.values(),
        #     ms1_res=M.ms1_resolution,
        #     ),
        # expand('process/{design}_iRT_{bin_iRT}_mz_{bin_mz}.barcode_ions.csv',
        #     design=M.name,
        #     bin_iRT=M.iRT_bin_names.values(),
        #     bin_mz=M.precursor_bin_names.values())
        # expand('process/{design}_iRT_{bin_iRT}_ms1_{ms1_range}.barcode_ions.csv',
        #     design=M.name,
        #     bin_iRT=list(M.iRT_bin_names.values()),
        #     ms1_range=list(M.ms1_selection_ranges.keys()))


rule generate_peptides:
    """Generate random peptides following composition rules. Save by precursor mz bin.
    """
    output:
        expand('process/{{design}}_{{run}}_{bin_mz}.peptides.csv', 
            bin_mz=M.precursor_bin_names.values())
    run:
        if M.rule_set == 'RJ_76':
            df_peptides = fly.design.permute_precursors(M.known_peptides, 
                M.num_permutations)
        elif 'exhaustive' in dir(M):
            prefix = wildcards['run']
            df_peptides = fly.design.generate_all_peptides(
                M.min_length, M.max_length, prefix, M.rule_set)
        else:
            df_peptides = fly.generate_peptide_set(
                M.num_to_generate, M.min_length, 
                M.max_length, M.rule_set, seed=int(wildcards.run))

        if 'peptide_gate' in dir(M):
            num_gen = len(df_peptides)
            df_peptides = fly.design.apply_peptide_gate(df_peptides, M.peptide_gate)
            msg = f'Gated {len(df_peptides):,}/{num_gen:,} generated peptides'
            print(msg, file=sys.stderr)

        df_peptides = (df_peptides
            .loc[lambda x: ~x['sequence'].str.contains(M.exclude_regex)]
            .assign(mz_bin=lambda x: 
                x['mz'].pipe(fly.bin_by_value, M.precursor_bins, 
                    M.precursor_bin_width))
            .query('mz_bin == mz_bin')
            .assign(mz_bin=lambda x: x['mz_bin'].map(M.precursor_bin_names))
            .assign(run=RUN_NAME)
        )

        # need to write empty csv files for empty bins
        for f, mz_bin in zip(output, M.precursor_bin_names.values()):
            (df_peptides.query('mz_bin == @mz_bin')
                .to_csv(f, index=None)
            )


rule predict_prosit:
    """Predict iRT and fragmentation spectra of precursor peptides.
    """
    input:
        expand('process/{{design}}_{run}_{{bin_mz}}.peptides.csv', 
            run=RUNS)
    output:
        expand('process/{{design}}_iRT_{bin_iRT}_mz_{{bin_mz}}.precursors.csv', 
            bin_iRT=M.iRT_bin_names.values())
    resources:
        gpu_mem_tenths=1
    run:
        d_spectra, d_irt = fly.load_prosit_models_cpu(
            MODEL_SPECTRA, MODEL_IRT)
        df_peptides = pd.concat([pd.read_csv(f) for f in input])
        if len(df_peptides) == 0:
            for f in output:
                pd.DataFrame().to_csv(f, index=None)
        else:
            df_predicted = (df_peptides
            # shuffle
            .sample(frac=1, replace=False, random_state=0)
            .head(M.pred_barcodes_per_mz_bin)
            .pipe(fly.add_prosit, d_spectra, d_irt, 
                M.normalized_collision_energy)
            .assign(iRT_bin=lambda x: 
                    x['iRT'].pipe(fly.bin_by_value, 
                                M.iRT_bins, M.iRT_bin_width))
            .pipe(fly.sort_by_spectral_efficiency)
            )
            df_predicted = df_predicted.query('iRT_bin == iRT_bin')
            # df_predicted.to_csv('predicted.csv', index=None)
            # print('iRT_bins', M.iRT_bin_names.values())


            for f, iRT_bin in zip(output, M.iRT_bin_names.values()):
                iRT_bin = float(iRT_bin)
                (df_predicted.query('iRT_bin == @iRT_bin')
                    .to_csv(f, index=None)
                    )


rule filter_barcodes:
    """Filter precursors to ensure each precursor generates ions that are unique
    within a precursor mz and iRT bin.
    """
    input:
        'process/{design}_iRT_{bin_iRT}_mz_{bin_mz}.precursors.csv'
    output:
        'process/{design}_iRT_{bin_iRT}_mz_{bin_mz}.barcode_ions.csv'        
    run:
        df_peptides = (pd.read_csv(input[0])
         .head(M.input_barcodes_per_iRT_mz_bin))

        if len(df_peptides) == 0:
            pd.DataFrame().to_csv(output[0], index=None)
        else:
            df_ions_selected, df_wide, barcode_ix = fly.snake_select_barcodes(
                df_peptides, M)
            df_ions_selected.to_csv(output[0], index=None)
            ax = fly.plot_ion_usage(df_wide, barcode_ix, M.ion_bins)
            ax.figure.savefig(output[0].replace('csv', 'png'), dpi=300)


rule filter_barcodes_ms1_range:
    """Original MS1 filtering strategy. Peptides are identified by iRT bin and precursor mz. 
    Precursors also generate unique ions within an mz range.
    """
    input:
        expand_ms1_range
    output:
        'process/{design}_iRT_{bin_iRT}_ms1_{ms1_range}.barcode_ions.csv'        
    run:
        df_peptides = (csv_frame(input, sort=False)
         .assign(iRT_dist=b'abs(iRT - iRT_bin)')
         .sort_values('iRT_dist')
         .drop('iRT_dist', axis=1)
         .groupby('mz_bin')
         .head(M.ms1_selection_input_max)
         .sort_values('mz_bin')
        )

        if len(df_peptides) == 0:
            pd.DataFrame().to_csv(output[0], index=None)
        else:
            df_ions_selected, df_wide, barcode_ix = fly.snake_select_barcodes(
                df_peptides, M, unique_col='mz_bin')
            (df_ions_selected
             .assign(ms1_range=wildcards.ms1_range)
             .to_csv(output[0], index=None)
            )


rule filter_by_ms1_resolution:
    """Simplified selection for high-resolution MS1 quantification. Require unique
    iRT and MS1 mz at specified resolution. Only allow barcodes with 
    `M.min_usable_ions` passing `M.usable_ion_gate`.
    """
    input:
        expand('process/{{design}}_iRT_{{bin_iRT}}_mz_{bin_mz}.precursors.csv', 
            bin_mz=M.precursor_bin_names.values())
    output:
        'process/{design}_iRT_{bin_iRT}.ms1_{ms1_res}.csv'
    run:
        df_barcodes = (csv_frame(input)
        .fillna(0) # some ions might be missing
        .pipe(fly.ms.filter_by_standards)
        )
        if len(df_barcodes) == 0:
            pd.DataFrame().to_csv(output[0], index=None)
        else:
            (df_barcodes
            .pipe(fly.design.add_usable_ion_count, M)
            .assign(usable_ion_count=lambda x: 
                    x['usable_ion_count'].clip(upper=M.min_usable_ions))
            .sort_values('usable_ion_count', ascending=False)
            .pipe(fly.design.add_mz_resolution_bins, int(wildcards['ms1_res']))
            .query('mz_res_bin_even')
            .drop_duplicates(['mz_res_bin_center', 'iRT_bin'])
            .to_csv(output[0], index=None)
            )


rule complete_ms1_resolution:
    """Consolidate 
    """
    input:
        expand('process/{design}_iRT_{bin_iRT}.ms1_{{ms1_res}}.csv',
            design=M.name,
            bin_iRT=M.iRT_bin_names.values())
    output:
        table='barcodes_ms1_{ms1_res}.csv',
        delta_mz='figures/ms1_{ms1_res}_delta_mz.png',
        iRT_vs_mz='figures/iRT_{ms1_res}_vs_mz.png',
    run:
        import seaborn as sns
        import matplotlib.pyplot as plt

        df_barcodes = csv_frame(input)
        cols = ['sequence', 'mz', 'mz_res_bin_center','iRT', 'iRT_bin', 'usable_ion_count']
        df_barcodes[cols].to_csv(output.table, index=None)

        os.makedirs('figures', exist_ok=True)
        nearest_mz = df_barcodes.set_index('sequence').sort_values('mz').groupby('iRT_bin')['mz'].diff()
        nearest_ratio = nearest_mz / list(df_barcodes['mz'])
        ax = nearest_ratio.dropna().pipe(np.log10).clip(upper=-4).hist()
        ax.set_xlabel('log10(delta mz / mz)')
        ax.figure.savefig(output.delta_mz)

        fg = (df_barcodes
        .assign(length=lambda x: x['sequence'].str.len())
        .pipe(sns.FacetGrid, height=8, aspect=2, hue='length')
        .map(plt.scatter, 'mz', 'iRT', s=4)
        .add_legend()
        .savefig(output.iRT_vs_mz)
        )

rule export_all_barcodes:
    """Generate a master barcode table
    """
    input:
        expand('process/{design}_iRT_{bin_iRT}_mz_{bin_mz}.precursors.csv', 
            design=M.name, 
            bin_mz=M.precursor_bin_names.values(),
            bin_iRT=M.iRT_bin_names.values(),
            )
    output:
        'all_barcodes.csv'
    run:
        cols = ['sequence', 'mz','iRT', 'usable_ion_count']
        df_barcodes = (csv_frame(input)
         .fillna(0) # some ions might be missing
         .pipe(fly.design.add_usable_ion_count, M)
         [cols].to_csv(output[0], index=None)
        )

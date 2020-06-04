sys.path.append(os.path.join(os.environ['HOME'], 'packages/prosit'))

GPU_MEM_FRACTION = 0.1

# IMPORTS

import postdoc.flycodes as fly
from postdoc.flycodes import designs
from postdoc.utils import timestamp

import pandas as pd
import inspect
METADATA = dict(inspect.getmembers(designs, inspect.isclass))[config['design']]

RUN_NAME = timestamp(METADATA.name)

MODEL_IRT = ('/home/dfeldman/flycodes/prosit_models/'
         'model_irt_prediction/')
MODEL_SPECTRA = ('/home/dfeldman/flycodes/prosit_models/'
                 'model_fragmentation_prediction/')

# RULES

RUNS = ['{:03d}'.format(x) for x in range(METADATA.num_generation_runs)]


rule all:
    input: 
        expand('process/{design}_{run}_{bin_mz}.peptides.csv', 
            design=METADATA.name,
            run=RUNS,
            bin_mz=METADATA.precursor_bin_names.values()),
        # expand('process/{design}_iRT_{bin_iRT}_mz_{bin_mz}.precursors.csv',
        #     design=METADATA.name,
        #     bin_iRT=METADATA.iRT_bin_names.values(),
        #     bin_mz=METADATA.precursor_bin_names.values())
        expand('process/{design}_iRT_{bin_iRT}_mz_{bin_mz}.barcode_ions.csv',
            design=METADATA.name,
            bin_iRT=METADATA.iRT_bin_names.values(),
            bin_mz=METADATA.precursor_bin_names.values())


rule generate_peptides:
    output:
        expand('process/{{design}}_{{run}}_{bin_mz}.peptides.csv', 
            bin_mz=METADATA.precursor_bin_names.values())
    run:
        seed = hash(wildcards['run']) % 2**32
        df_precursors_all = (fly.generate_peptide_set(
            METADATA.num_to_generate, METADATA.min_length, METADATA.max_length,
            METADATA.rule_set, seed=seed)
            .assign(mz_bin=lambda x: 
                x['mz'].pipe(fly.bin_by_value, METADATA.precursor_bins, 
                    METADATA.precursor_bin_width))
             .query('mz_bin == mz_bin')
             .assign(mz_bin=lambda x: x['mz_bin'].map(METADATA.precursor_bin_names))
             .assign(run=RUN_NAME)
        )

        # need to write empty csv files for empty bins
        for f, mz_bin in zip(output, METADATA.precursor_bin_names.values()):
            (df_precursors_all.query('mz_bin == @mz_bin')
                .to_csv(f, index=None)
            )


rule predict_prosit:
    input:
        expand('process/{{design}}_{run}_{{bin_mz}}.peptides.csv', 
            run=RUNS)
    output:
        expand('process/{{design}}_iRT_{bin_iRT}_mz_{{bin_mz}}.precursors.csv', 
            bin_iRT=METADATA.iRT_bin_names.values())
    resources:
        gpu_mem_tenths=1
    run:
        d_spectra, d_irt = fly.load_prosit_models(
            MODEL_IRT, MODEL_SPECTRA, gpu_mem_fraction=GPU_MEM_FRACTION)
        df_precursors = (pd.concat([pd.read_csv(f) for f in input])
            .rename(columns={'orig_seq': 'sequence'})
            )

        num_to_predict = min(
            len(df_precursors), 
            METADATA.pred_barcodes_per_mz_bin)
        
        df_predicted = (df_precursors
         .sample(frac=1, replace=False, random_state=0)
         .head(num_to_predict)
         .pipe(fly.add_prosit, d_spectra, d_irt, 
            METADATA.normalized_collision_energy)
         .assign(iRT_bin=lambda x: 
                 x['iRT'].pipe(fly.bin_by_value, 
                               METADATA.iRT_bins, METADATA.iRT_bin_width))
         .query('iRT_bin == iRT_bin')
         .pipe(fly.sort_by_spectral_efficiency)
        )

        for f, iRT_bin in zip(output, METADATA.iRT_bin_names.values()):
            (df_predicted.query('iRT_bin == @iRT_bin')
                .head(METADATA.input_barcodes_per_iRT_mz_bin)
                .to_csv(f, index=None)
                )


rule filter_barcodes:
    input:
        'process/{design}_iRT_{bin_iRT}_mz_{bin_mz}.precursors.csv'
    output:
        'process/{design}_iRT_{bin_iRT}_mz_{bin_mz}.barcode_ions.csv'        
    run:
        df_peptides = pd.read_csv(input[0])

        if len(df_peptides) == 0:
            pd.DataFrame().to_csv(output[0], index=None)
        else:
            df_ions_selected, df_wide, barcode_ix = fly.snake_select_barcodes(
                df_peptides, METADATA)
            df_ions_selected.to_csv(output[0], index=None)
            ax = fly.plot_ion_usage(df_wide, barcode_ix, METADATA.ion_bins)
            ax.figure.savefig(output[0].replace('csv', 'png'), dpi=300)


"""
squeue --user=dfeldman

sbatch -p gpu --mem=80g --gres=gpu:rtx2080:1 -c 10 run_003.sh
sbatch -p gpu --mem=80g --gres=gpu:rtx2080:1 -c 10 run_004.sh

conda activate prosit5
cd /home/dfeldman/flycodes/run_003
snakemake -k -s /home/dfeldman/packages/postdoc/scripts/flycodes.smk \
 --config design=DESIGN_3 --resources gpu_mem_tenths=6 --cores

conda activate prosit5
cd /home/dfeldman/flycodes/run_003
snakemake -k -s /home/dfeldman/packages/postdoc/scripts/flycodes.smk \
 --config design=DESIGN_4 --resources gpu_mem_tenths=6 --cores

"""
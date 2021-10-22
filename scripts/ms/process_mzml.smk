from postdoc.flycodes import ms_app
import pandas as pd
import os

config = ms_app.load_config()['process_mzml']

df_samples = pd.read_csv(ms_app.sample_table)
SAMPLES = list(df_samples['sample'])
os.environ['PATH'] += ':' + ms_app.openms_dir
os.environ['PATH'] += ':' + ms_app.tpp_dir



rule all:
    input:
        # expand('{sample}.calibratedRT.mzML', sample=SAMPLES),
        # expand('{sample}.calibratedRT.pepXML', sample=SAMPLES),
        # expand('{sample}.calibrated.pepXML', sample=SAMPLES),
        # expand('{sample}.calibrated.mzML', sample=SAMPLES),
        # expand('{sample}.filtered.pepXML', sample=SAMPLES),
        expand(f'{{sample}}.{config["request"]}.mzML', sample=SAMPLES),
        expand(f'{{sample}}.{config["request"]}.pepXML', sample=SAMPLES),
        expand(f'ms1/{{sample}}.{config["request"]}.ms1.hdf', sample=SAMPLES),

rule limit_range:
    input: 
        '{sample}.mzML'
    output: 
        '{sample}.filtered.mzML',
    params:
        flags=config['mzML_filter'],
        app=ms_app.ms_app,
    shell:
        """
        FileFilter -in {input} -out {output[0]} {params.flags}
        """

rule comet:
    input:
        '{sample}.{tag}.mzML'
    output:
        '{sample}.{tag}.pepXML',
        '{sample}.{tag}.idXML',
    params:
        database='../' + ms_app.barcodes_by_design,
        comet=config['comet_params'],
        idXML_filter=config['idXML_filter'],
    shell:
        """
        comet -D{params.database} -P{params.comet} {input}
        mv {wildcards.sample}.{wildcards.tag}.pep.xml {output[0]}
        IDFileConverter -in {output[0]} -out {output[1]}
        IDFilter -in {output[1]} -out {output[1]} {params.idXML_filter}
        """

rule calibrate_mz:
    input:
        '{sample}.filtered.mzML',
        '{sample}.filtered.idXML',
    output:
        '{sample}.calibrated.mzML'
    params:
        flags=config['mz_calibration']
    shell:
        """ExternalCalibration -in {input[0]} -out {output} {params.flags}
        """

rule calibrate_rt:
    input:
        expand('{sample}.calibrated.idXML', sample=SAMPLES),
    output:
        expand('{sample}.calibratedRT.trafoXML', sample=SAMPLES),
    params:
        flags=f'{config["rt_calibration"]} -reference:file {config["rt_reference"]}.calibrated.idXML'
    shell:
        """
        MapAlignerIdentification -in {input} -trafo_out {output} {params.flags}
        """
   
rule apply_rt_calibration:
    input:
        '{sample}.calibrated.mzML',
        '{sample}.calibratedRT.trafoXML',
    output:
        '{sample}.calibratedRT.mzML'
    shell:
        """
        MapRTTransformer -in {input[0]} -trafo_in {input[1]} -out {output}
        """

rule export_ms1_filtered:
    input:
        '{sample}.filtered.mzML'
    output:
        
    run:
        os.makedirs('ms1', exist_ok=True)
        ms_app.export_ms1(input[0]).to_csv(output[0], index=None)


rule export_ms1:
    input:
        '{sample}.{tag}.mzML'
    output:
        'ms1/{sample}.{tag}.ms1.hdf'
    run:
        os.makedirs('ms1', exist_ok=True)
        ms_app.export_ms1(input[0]).to_hdf(output[0], 'x', mode='w')

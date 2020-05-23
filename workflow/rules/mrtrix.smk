# general rule file for mrtrix processing
# requires prepdwi and fmriprep to be already run

rule convert_to_mif:
    input:
        dwi_nii = join(config['prepdwi_dir'],config['in_dwi_nii']),
        dwi_bval = join(config['prepdwi_dir'],config['in_dwi_bval']),
        dwi_bvec = join(config['prepdwi_dir'],config['in_dwi_bvec']),
        mask_nii = join(config['prepdwi_dir'],config['in_dwi_mask_nii'])
    output:
        #use simplified naming for mrtrix work files
        dwi_mif = 'results/sub-{subject}/mrtrix/dwi.mif',
        b0_mif = 'results/sub-{subject}/mrtrix/b0.mif',
        mask_mif = 'results/sub-{subject}/mrtrix/mask.mif'
    log:
        'results/logs/convert_to_mif/sub-{subject}.log'
    envmodules: 
        'mrtrix'
    #group: 'tracking'
    container: config['singularity']
    shell:
        "(mrconvert {input.dwi_nii} {output.dwi_mif} -fslgrad {input.dwi_bvec} {input.dwi_bval} -datatype float32 -strides 0,0,0,1 &&"
        "mrconvert {input.mask_nii} {output.mask_mif} &&"
        "dwiextract {output.dwi_mif} - -bzero | mrmath - mean {output.b0_mif} -axis 3) &> {log}"


rule gen_5tt:
    input: 
        lut_input = 'resources/lut/FreeSurferColorLUT.txt',
        lut_output = 'resources/lut/FreeSurfer2ACT.txt',
        aparcaseg_nii = config['in_aparcaseg_nii']
    output:
        indices_mif = 'results/sub-{subject}/mrtrix/indices.mif',
        cgm_mif = temp('results/sub-{subject}/mrtrix/cgm.mif'),
        sgm_mif = temp('results/sub-{subject}/mrtrix/sgm.mif'),
        wm_mif = temp('results/sub-{subject}/mrtrix/wm.mif'),
        csf_mif = temp('results/sub-{subject}/mrtrix/csf.mif'),
        path_mif = temp('results/sub-{subject}/mrtrix/path.mif'),
        seg_5tt = 'results/sub-{subject}/mrtrix/seg_5tt.mif'
    log: 'results/logs/gen_5tt/sub-{subject}.log'
    #group: 'tracking'
    container: config['singularity']
    shell:
        "(labelconvert {input.aparcaseg_nii} {input.lut_input} {input.lut_output} {output.indices_mif} && "
        "mrcalc {output.indices_mif} 1 -eq {output.cgm_mif} &&"
        "mrcalc {output.indices_mif} 2 -eq {output.sgm_mif} &&"
        "mrcalc {output.indices_mif} 3 -eq {output.wm_mif} &&"
        "mrcalc {output.indices_mif} 4 -eq {output.csf_mif} &&"
        "mrcalc {output.indices_mif} 5 -eq {output.path_mif} &&"
        "mrcat {output.cgm_mif} {output.sgm_mif} {output.wm_mif} {output.csf_mif} {output.path_mif} - -axis 3 | mrconvert - {output.seg_5tt} -datatype float32 && "
        "5ttcheck {output.seg_5tt}) &> {log}"
"""

rule gen_5tt_fsl:
    input:
        t1 = join(config['fmriprep_dir'],'sub-{subject}/anat/sub-{subject}_desc-preproc_T1w.nii.gz'),
        mask = join(config['fmriprep_dir'],'sub-{subject}/anat/sub-{subject}_desc-brain_mask.nii.gz')
    output:
        seg_5tt = 'results/sub-{subject}/mrtrix/seg_5tt.mif',
        tempdir = directory('results/mrtrix/sub-{subject}/5ttgen')
    threads: 8
    log: 'results/logs/gen_5tt/sub-{subject}.log'
    envmodules: 
        'mrtrix',
        'fsl'
    #group: 'tracking'
    shell:
        '5ttgen fsl {input.t1} {output.seg_5tt} -mask {input.mask} -nocrop -tempdir {output.tempdir} -nocleanup -nthreads {threads} &> {log} && ' 
        '5ttcheck {output.seg_5tt} &> {log}'
"""
rule estimate_response:
    input:         
        dwi_mif = 'results/sub-{subject}/mrtrix/dwi.mif',
        seg_5tt = 'results/sub-{subject}/mrtrix/seg_5tt.mif',
    params:
        tempdir = 'results/sub-{subject}/tmp_dwi2response'
    #can add additional params here (lmax)
    output:
        rf_wm = 'results/sub-{subject}/mrtrix/rf_wm.txt',
        rf_gm = 'results/sub-{subject}/mrtrix/rf_gm.txt',
        rf_csf = 'results/sub-{subject}/mrtrix/rf_csf.txt',
        rf_voxels = 'results/sub-{subject}/mrtrix/rf_voxels.mif'
    threads: 8
    envmodules: 
        'mrtrix'
    container: config['singularity']
    log:
        'results/logs/estimate_response/sub-{subject}.log'

    shell:
        "dwi2response -tempdir {params.tempdir} msmt_5tt {input.dwi_mif} {input.seg_5tt} {output.rf_wm} {output.rf_gm} {output.rf_csf} -voxels {output.rf_voxels} -nthreads {threads} &> {log}"
    
rule compute_fod:
    input:
        dwi_mif = 'results/sub-{subject}/mrtrix/dwi.mif',
        rf_wm = 'results/sub-{subject}/mrtrix/rf_wm.txt',
        rf_gm = 'results/sub-{subject}/mrtrix/rf_gm.txt',
        rf_csf = 'results/sub-{subject}/mrtrix/rf_csf.txt',
        mask_mif = 'results/sub-{subject}/mrtrix/mask.mif'
    output:
        fod_wm = 'results/sub-{subject}/mrtrix/fod_wm.mif',
        fod_gm = 'results/sub-{subject}/mrtrix/fod_gm.mif',
        fod_csf = 'results/sub-{subject}/mrtrix/fod_csf.mif'
    resources:
        time = 6*60, #in minutes
    threads: 8
    envmodules: 
        'mrtrix'
    container: config['singularity']
    log:
        'results/logs/compute_fod/sub-{subject}.log'
    shell:
        "dwi2fod -mask {input.mask_mif} msmt_csd {input.dwi_mif} {input.rf_wm}  {output.fod_wm} {input.rf_gm}  {output.fod_gm} {input.rf_csf}  {output.fod_csf} &> {log}"

rule mtnormalise:
    input:
        fod_wm = 'results/sub-{subject}/mrtrix/fod_wm.mif',
        fod_gm = 'results/sub-{subject}/mrtrix/fod_gm.mif',
        fod_csf = 'results/sub-{subject}/mrtrix/fod_csf.mif',
        mask_mif = 'results/sub-{subject}/mrtrix/mask.mif'
    output:
        fodn_wm = 'results/sub-{subject}/mrtrix/fodn_wm.mif',
        fodn_gm = 'results/sub-{subject}/mrtrix/fodn_gm.mif',
        fodn_csf = 'results/sub-{subject}/mrtrix/fodn_csf.mif'
    envmodules: 
        'mrtrix'
    container: config['singularity']
    log:
        'results/logs/mtnormalise/sub-{subject}.log'
    #group: 'tracking'
    shell:
        "mtnormalise -mask {input.mask_mif} {input.fod_wm} {output.fodn_wm}  {input.fod_gm} {output.fodn_gm}  {input.fod_csf} {output.fodn_csf}  &> {log}"

rule gen_tracts:
    input:
        fodn_wm = 'results/sub-{subject}/mrtrix/fodn_wm.mif',
    output:
        tracts = 'results/sub-{subject}/mrtrix/tracts.tck'
    params:
        ntracts = f"{config['ntracts_million']}M"
    resources:
        time = int(config['ntracts_million'] * 0.25 * 60),
        mem_mb = 4000 #in mb
    threads: 8 # num cores
    envmodules: 
        'mrtrix'
    container: config['singularity']
    log:
        'results/logs/gen_tracts/sub-{subject}.log'
    #group: 'tracking'
    shell:
        "tckgen {input.fodn_wm} {output.tracts} -select {params.ntracts} -seed_dynamic {input.fodn_wm} -backtrack -crop_at_gmwmi -maxlength 250 -cutoff 0.06 -nthreads {threads} &> {log}"
   

rule run_sift2:
    input:
        tracts = 'results/sub-{subject}/mrtrix/tracts.tck',
        fodn_wm = 'results/sub-{subject}/mrtrix/fodn_wm.mif',
        seg_5tt = 'results/sub-{subject}/mrtrix/seg_5tt.mif'
    output:
        weights = 'results/sub-{subject}/mrtrix/sift2_weights.txt'
    threads: 8 # num cores
    envmodules: 
        'mrtrix'
    container: config['singularity']
    log:
        'results/logs/run_sift2/sub-{subject}.log'
    #group: 'tracking'
    shell:
        "tcksift2 -act {input.seg_5tt} {input.tracts} {input.fodn_wm} {output.weights} -nthreads {threads} &> {log}"

rule convert_atlas_labels:
    input:
        aparcaseg_nii = config['in_aparcaseg_nii'],
        lut_in = 'resources/lut/FreeSurferColorLUT.txt',
        lut_out = 'resources/lut/fs_default.txt'
    output:
        atlas_labels = 'results/sub-{subject}/mrtrix/atlas_aparcaseg.nii.gz'
    envmodules: 
        'mrtrix'
    container: config['singularity']
    log:
        'results/logs/convert_atlas_labels/sub-{subject}.log'
    #group: 'tracking'
    shell:
        "labelconvert {input.aparcaseg_nii} {input.lut_in} {input.lut_out} {output.atlas_labels} &> {log}"

rule gen_connectome:
    input:  
        tracts = 'results/sub-{subject}/mrtrix/tracts.tck',
        atlas_labels = 'results/sub-{subject}/mrtrix/atlas_aparcaseg.nii.gz',
        weights = 'results/sub-{subject}/mrtrix/sift2_weights.txt'
    output:
        connectome = 'results/sub-{subject}/mrtrix/connectome_aparcaseg.csv'
    log:
        'results/logs/gen_connectome/sub-{subject}.log'
    container: config['singularity']
    shell:
        "tck2connectome -tck_weights_in {input.weights} {input.tracts} {input.atlas_labels} {output.connectome} &> {log}"

rule visualize_connectome:
    input:
        connectome = 'results/sub-{subject}/mrtrix/connectome_aparcaseg.csv'
    output:
        report(join('plots','sub-{subject}_connectome_vis.png'),caption='../report/connectome_vis.rst',category='Connectome Visualization')
    conda: '../envs/vis.yml'
    script: '../scripts/visualize_connectome.py'


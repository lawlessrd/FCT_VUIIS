#!/usr/bin/env bash
#
# This script allows testing the compiled matlab, assuming the correct matlab
# runtime is installed. Much better to make sure it's working before actually 
# building the singularity container.

# Single fmri input
bin/run_spm12.sh /usr/local/MATLAB/MATLAB_Runtime/v97 \
    function fct_entrypoint \
    t1_niigz "/home/dylan/Documents/FCT/INPUTS/T1w.nii.gz" \
    fmri_niigz "/home/dylan/Documents/FCT/INPUTS/fmri.nii.gz" \
    out_dir "/home/dylan/Documents/FCT/OUTPUTS" \
    xnat_project "OASIS-3" \
    xnat_subject "TEST_SUBJ" \
    xnat_session "TEST_SESS"
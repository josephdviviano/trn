## Dual Phase Encoded Stimuli Mean Generation w. AFNI & MATLAB
## Allows for multi-session registrations

## June 19th 2012, Joseph Viviano
## 
## SESS1 ... SESSn folders, func files 
##
# rename functional files in /SESSn

# This script:
# 1) Despike
# 2) Linear detrend
# 3) Records timeseries mean
# 4) Convert to % signal change
# 5) Record nusiance time series
# 6) 3dDetrend:
#    a) detrend (polort = 2)
#    b) regress nusiances
# 7) Smooth within mask

# Output: filterAll(global, WM, ventricles, motion) + smoothed
#         filterVent(ventricles, motion) + smoothed
#         filterHead(motion) + smoothed
#         filterNone() + smoothed

## SETTINGS ###################################################################
SUBJECTS="1 2 5 7 9 11"            # list of subject numbers to analyze
DIREXP='flicker_120_thal'

CORES=8         # number of jobs (should equal the number of cores you have)
FWHM=2          # spatial smoothing (mm)
###############################################################################
for SUB in $SUBJECTS; do

SESS=`ls -l -d ${DATA_DIR}/s${SUB}/${DIREXP}/SESS* | wc -l` 
for ((s=1;s<=SESS;s++)); do

## set number of runs for current session
RUN=`ls -l ${DATA_DIR}/s${SUB}/${DIREXP}/SESS$s/raw/*func* | wc -l`
rm /tmp/list.txt
for ((i=1;i<=RUN;i++)); do echo $i >> /tmp/list.txt; done

    # despike
    if [ ! -f ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_despike1.nii.gz ]; then
    cat /tmp/list.txt | parallel -P $CORES \
    3dDespike -prefix ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_despike{}.nii.gz \
              ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_motion{}.nii.gz
    fi
     
    # Convert to % signal change
    if [ ! -f ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_pct1.nii.gz ]; then
    cat /tmp/list.txt | parallel -P $CORES \
    3dDetrend -prefix ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_tmp_det{}.nii.gz \
              -polort 1 ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_despike{}.nii.gz

    cat /tmp/list.txt | parallel -P $CORES \
    3dTstat -prefix ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_tmp_mean{}.nii.gz \
            -mean ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_despike{}.nii.gz

    cat /tmp/list.txt | parallel -P $CORES \
    3dTstat -prefix ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_tmp_detMean{}.nii.gz \
            -mean ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_tmp_det{}.nii.gz

    cat /tmp/list.txt | parallel -P $CORES \
    3dCalc -prefix ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_tmp_detrend{}.nii.gz \
           -a ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_tmp_det{}.nii.gz \
           -b ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_tmp_mean{}.nii.gz \
           -c ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_tmp_detMean{}.nii.gz \
           -expr "'a+b-c'"

    cat /tmp/list.txt | parallel -P $CORES \
    3dTstat -prefix ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_tmp_pctMean{}.nii.gz \
            -mean ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_tmp_detrend{}.nii.gz

    cat /tmp/list.txt | parallel -P $CORES \
    3dcalc -prefix ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_pct{}.nii.gz \
           -a ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_tmp_detrend{}.nii.gz \
           -b ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_tmp_pctMean{}.nii.gz \
           -expr "'(a-b)/b * 100'"

    purge
    rm ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_tmp*
    fi

    ## Compute all nusiance time series
    for ((i=1;i<=RUN;i++)); do 
    
        # Global mean
        if [ ! -f ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/params_globalMean${i}.1D ]; then
        3dmaskave -mask ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/anat_EPI_mask_reg.nii.gz \
                  -quiet ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_pct${i}.nii.gz \
                  > ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/params_globalMean${i}.1D
        fi

        # Global mean derivative
        if [ ! -f ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/params_globalDeriv${i}.1D ]; then
        1d_tool.py -infile ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/params_globalMean${i}.1D \
                   -derivative -overwrite \
                   -write ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/params_globalDeriv${i}.1D
        fi

        # Ventricles
        if [ ! -f ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/params_ventricleMean${i}.1D ]; then
        3dmaskave -mask ${DATA_DIR}/s${SUB}/${DIREXP}/anat_vent_mask.nii.gz \
                  -quiet ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_pct${i}.nii.gz \
                  > ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/params_ventricleMean${i}.1D
        fi

        # White matter
        if [ ! -f ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/params_WMMean${i}.1D ]; then
        3dmaskave -mask ${DATA_DIR}/s${SUB}/${DIREXP}/anat_WM_mask.nii.gz \
                  -quiet ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_pct${i}.nii.gz \
                  > ${DATA_DIR}/s${SUB}/${DIREXP}/SESS$s/params_WMMean${i}.1D
        fi

    done

    # Detrend (Polort = 2), regress nusiances [full, ventricle, none], smooth
    if [ ! -f ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_filterFull1.nii.gz ]; then
    cat /tmp/list.txt | parallel -P $CORES \
    3dDetrend -prefix ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_filterFull{}.nii.gz \
              -polort 2 \
              -vector ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/params_ventricleMean{}.1D \
              -vector ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/params_mot{}.1D \
              -vector ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/params_globalMean{}.1D \
              -vector ${DATA_DIR}/s${SUB}/${DIREXP}/SESS$s/params_WMMean{}.1D \
              ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_pct{}.nii.gz
    fi
    
    if [ ! -f ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_filterVent1.nii.gz ]; then
    cat /tmp/list.txt | parallel -P $CORES \
    3dDetrend -prefix ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_filterVent{}.nii.gz \
              -polort 2 \
              -vector ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/params_ventricleMean{}.1D \
              -vector ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/params_mot{}.1D \
              -vector ${DATA_DIR}/s${SUB}/${DIREXP}/SESS$s/params_WMMean{}.1D \
              ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_pct{}.nii.gz
    fi

    if [ ! -f ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_filterHead1.nii.gz ]; then
    cat /tmp/list.txt | parallel -P $CORES \
    3dDetrend -prefix ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_filterHead{}.nii.gz \
              -polort 2 \
              -vector ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/params_mot{}.1D \
              ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_pct{}.nii.gz
    fi

    if [ ! -f ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_filterNone1.nii.gz ]; then
    cat /tmp/list.txt | parallel -P $CORES \
    3dDetrend -prefix ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_filterNone{}.nii.gz \
              -polort 2 \
              ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_pct{}.nii.gz
    fi

    if [ ! -f ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_smoothFull1.nii.gz ]; then
    cat /tmp/list.txt | parallel -P $CORES \
    3dBlurInMask -prefix ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_smoothFull{}.nii.gz \
                 -mask ${DATA_DIR}/s${SUB}/${DIREXP}/anat_THAL_mask.nii.gz \
                 -FWHM ${FWHM} \
                 -input ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_filterFull{}.nii.gz
    fi        

    if [ ! -f ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_SmoothVent1.nii.gz ]; then
    cat /tmp/list.txt | parallel -P $CORES \
    3dBlurInMask -prefix ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_smoothVent{}.nii.gz \
                 -mask ${DATA_DIR}/s${SUB}/${DIREXP}/anat_THAL_mask.nii.gz \
                 -FWHM ${FWHM} \
                 -input ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_filterVent{}.nii.gz
    fi

    if [ ! -f ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_smoothHead1.nii.gz ]; then
     cat /tmp/list.txt | parallel -P $CORES \
     3dBlurInMask -prefix ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_smoothHead{}.nii.gz \
                  -mask ${DATA_DIR}/s${SUB}/${DIREXP}/mask_LGN_TRN_nopul.nii.gz \
                  -FWHM ${FWHM} \
                  -input ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_filterHead{}.nii.gz
    fi
    
    if [ ! -f ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_SmoothNone1.nii.gz ]; then
    cat /tmp/list.txt | parallel -P $CORES \
    3dBlurInMask -prefix ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_smoothNone{}.nii.gz \
                 -mask ${DATA_DIR}/s${SUB}/${DIREXP}/anat_THAL_mask.nii.gz \
                 -FWHM ${FWHM} \
                 -input ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_filterNone{}.nii.gz
    fi
    done
    purge
done

## JDV Jul 17 2013 ################################################################
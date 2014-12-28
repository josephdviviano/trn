## Dual Phase Encoded Stimuli Mean Generation w. AFNI & MATLAB
## Allows for multi-session registrations 
## SESS1 ... SESSn folders, with anat_WHEPI, anat_T1, anat_PD in master folder

# This script:
# 1) Delete initial time points
# 2) Time shifts data
# 3) Registrations: EPI --> WHEPI --> T1
 
## WARNING I CHANGED MOTION CORRECTION PARAMATERS 6 --> 12 --> 20 --> 50 might take forever ##
##                                        TWOBLUR 3 --> 2

## SETTINGS ###################################################################
CORES=8         # number of jobs (should equal the number of cores you have)
SUBJECTS="5"           # list of subject numbers to analyze
DIREXP='flicker_120_thal'

RES=0.75        # upsampled EPI resolution (iso)
REWINDTRS=8     # number of TRs to end of rewind
TR='1250ms'
## for each subject, set number of sessions and runs###########################

for SUB in $SUBJECTS; do

SESS=`ls -l -d ${DATA_DIR}/s${SUB}/${DIREXP}/SESS* | wc -l` 
for ((s=1;s<=SESS;s++)); do

RUN=`ls -l ${DATA_DIR}/s${SUB}/${DIREXP}/SESS$s/raw/*func* | wc -l`
rm /tmp/list.txt
for ((i=1;i<=RUN;i++)); do echo $i >> /tmp/list.txt; done

    ## delete rewind TRs
    if [ ! -f ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_tshift1.nii.gz ]; then
    cat /tmp/list.txt | parallel -P $CORES \
    3dcalc -prefix ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_tmp_del{}.nii.gz \
           -a ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/raw/func{}.nii.gz[${REWINDTRS}..$] \
           -expr 'a'

    ## slice time correction
    cat /tmp/list.txt | parallel -P $CORES \
    3dTshift -prefix ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_tshift{}.nii.gz -verbose -Fourier \
             -tpattern altplus -TR $TR ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_tmp_del{}.nii.gz
    fi

    ## deoblique all runs
    for ((i=1;i<=RUN;i++)); do
        if [ ! -f ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/mat_EPI_deoblique_${i}.aff12.1D ]; then
        3dWarp -deoblique -quintic -verb -prefix ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_tmp_ob${i}.nii.gz ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_tshift${i}.nii.gz > ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/mat_EPI_tmp_deoblique_${i}.1D
        sed '1 d' ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/mat_EPI_tmp_deoblique_${i}.1D > ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/mat_EPI_tmp_deoblique_headerstrip_${i}.1D
        cat_matvec -ONELINE ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/mat_EPI_tmp_deoblique_headerstrip_${i}.1D > ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/mat_EPI_deoblique_${i}.aff12.1D
        cat_matvec -ONELINE ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/mat_EPI_deoblique_${i}.aff12.1D -I > ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/mat_EPI_deoblique_inv_${i}.aff12.1D
        fi
    done
    rm ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/*mat_EPI_tmp_deoblique*

    ## compute between-run motion paramaters for registration & motion censoring
    cat /tmp/list.txt | parallel -P $CORES \
    3dVolreg -prefix ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_tmp_mot{}.nii.gz \
             -base ${DATA_DIR}'/s'${SUB}'/'${DIREXP}'/SESS'${s}'/func_tmp_ob1.nii.gz[8]' \
             -twopass -twoblur 2 -twodup -coarse 10 3 \
             -Fourier -zpad 50 -float \
             -1Dfile ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/params_mot{}.1D \
             -1Dmatrix_save ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/params_reg{}.aff12.1D \
              ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_tmp_ob{}.nii.gz

    ## create session EPI
    3dMean -prefix ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/anat_EPI_tmp_mean.nii.gz ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/*func_tmp_mot*
    3dTstat -prefix ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/anat_EPI_tmp_vol.nii.gz ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/anat_EPI_tmp_mean.nii.gz
    3dAutomask -prefix ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/anat_EPI_mask_${s}.nii.gz -clfrac 0.3 -SI 100 -peels 4 ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/anat_EPI_tmp_vol.nii.gz
    3dCalc -prefix ${DATA_DIR}/s${SUB}/${DIREXP}/anat_EPI_brain_${s}.nii.gz -a ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/anat_EPI_tmp_vol.nii.gz -b ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/anat_EPI_mask_${s}.nii.gz -expr 'a*b'
    rm ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/*EPI_tmp*

    ## create session Whole Head EPI
    3dWarp -deoblique -quintic -prefix ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/anat_WHEPI_tmp_ob.nii.gz ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/raw/anat_WHEPI.nii.gz
    3dResample -dxyz ${RES} ${RES} ${RES} \
               -rmode Cu \
               -prefix ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/anat_WHEPI_tmp_resize.nii.gz \
               -inset ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/anat_WHEPI_tmp_ob.nii.gz
    3dTstat -prefix ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/anat_WHEPI_tmp_vol.nii.gz ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/anat_WHEPI_tmp_resize.nii.gz
    bet ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/anat_WHEPI_tmp_vol.nii.gz ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/anat_WHEPI_tmp.nii.gz -R -m -n -f 0.3
    3dAutomask -prefix ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/anat_WHEPI_tmp_cropmask.nii.gz -clfrac 0.3 -SI 110 ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/anat_WHEPI_tmp_vol.nii.gz
    3dCalc -prefix ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/anat_WHEPI_brain.nii.gz \
           -a ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/anat_WHEPI_tmp_vol.nii.gz \
           -b ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/anat_WHEPI_tmp_mask.nii.gz \
           -c ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/anat_WHEPI_tmp_cropmask.nii.gz \
           -expr 'a*b*c'
    mv ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/anat_WHEPI_brain.nii.gz ${DATA_DIR}/s${SUB}/${DIREXP}/anat_WHEPI_brain_${s}.nii.gz
    rm ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/*WHEPI_tmp*

done

## register each session EPI to the session WHEPI
for ((s=1;s<=SESS;s++)); do
    cd ${DATA_DIR}/s${SUB}/${DIREXP}
    align_epi_anat.py -anat anat_WHEPI_brain_${s}.nii.gz \
                      -epi anat_EPI_brain_${s}.nii.gz \
                      -epi_base 0 -epi2anat -suffix EPI_to_WHEPI \
                      -anat_has_skull no -epi_strip None \
                      -volreg off -tshift off -deoblique off \
                      -giant_move -partial_coronal \
                      -cost lpa

    mv ${DATA_DIR}/s${SUB}/${DIREXP}/anat_EPI_brain_${s}.nii.gzEPI_to_WHEPI_mat.aff12.1D ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/mat_EPI_to_WHEPI_${s}.aff12.1D
    mv ${DATA_DIR}/s${SUB}/${DIREXP}/anat_WHEPI_brain_${s}.nii.gzEPI_to_WHEPI_mat.aff12.1D ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/mat_WHEPI_to_EPI_${s}.aff12.1D
    mv ${DATA_DIR}/s${SUB}/${DIREXP}/anat_EPI_brain_${s}.nii.gzEPI_to_WHEPI.nii.gz ${DATA_DIR}/s${SUB}/${DIREXP}/reg_EPI_to_WHEPI_${s}.nii.gz
done

## register each session WHEPI to T1
for ((s=1;s<=SESS;s++)); do
    cd ${DATA_DIR}/s${SUB}/${DIREXP}
    align_epi_anat.py -anat anat_T1_brain.nii.gz \
                      -epi anat_WHEPI_brain_${s}.nii.gz \
                      -epi_base 0 -epi2anat -suffix WHEPI_to_T1 \
                      -anat_has_skull no -epi_strip None \
                      -volreg off -tshift off -deoblique off \
                      -giant_move -cost lpc+zz

    mv ${DATA_DIR}/s${SUB}/${DIREXP}/anat_WHEPI_brain_${s}.nii.gzWHEPI_to_T1_mat.aff12.1D ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/mat_WHEPI_to_T1_${s}.aff12.1D
    mv ${DATA_DIR}/s${SUB}/${DIREXP}/anat_T1_brain.nii.gzWHEPI_to_T1_mat.aff12.1D ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/mat_T1_to_WHEPI_${s}.aff12.1D
    mv ${DATA_DIR}/s${SUB}/${DIREXP}/anat_WHEPI_brain_${s}.nii.gzWHEPI_to_T1.nii.gz ${DATA_DIR}/s${SUB}/${DIREXP}/reg_WHEPI_to_T1_${s}.nii.gz 
done

## create master EPI target
cat_matvec ${DATA_DIR}/s${SUB}/${DIREXP}/SESS1/mat_EPI_to_WHEPI_1.aff12.1D ${DATA_DIR}/s${SUB}/${DIREXP}/SESS1/mat_WHEPI_to_T1_1.aff12.1D > ${DATA_DIR}/s${SUB}/${DIREXP}/mat_regTarget.aff12.1D
3dWarp -matvec_out2in ${DATA_DIR}/s${SUB}/${DIREXP}/mat_regTarget.aff12.1D -prefix ${DATA_DIR}/s${SUB}/${DIREXP}/anat_EPI_tmp_regTarget.nii.gz -cubic -newgrid ${RES} -zpad 0 ${DATA_DIR}/s${SUB}/${DIREXP}/anat_EPI_brain_1.nii.gz
3dAutobox -prefix ${DATA_DIR}/s${SUB}/${DIREXP}/anat_EPI_regTarget.nii.gz -input ${DATA_DIR}/s${SUB}/${DIREXP}/anat_EPI_tmp_regTarget.nii.gz
rm ${DATA_DIR}/s${SUB}/${DIREXP}/*EPI_tmp*

## compute all funcs at 2x resolution in register with the T1
for ((s=1;s<=SESS;s++)); do
    mkdir ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/params 
    RUN=`ls -l ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/raw/func* | wc -l`
    for ((i=1;i<=RUN;i++)); do
        cat_matvec ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/mat_EPI_deoblique_inv_${i}.aff12.1D ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/params_reg${i}.aff12.1D ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/mat_EPI_to_WHEPI_${s}.aff12.1D ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/mat_WHEPI_to_T1_${s}.aff12.1D > ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/params/reg${i}.aff12.1D
        3dAllineate -prefix ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_motion${i}.nii.gz \
                    -input ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_tshift${i}.nii.gz \
                    -1Dmatrix_apply ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/params/reg${i}.aff12.1D \
                    -master ${DATA_DIR}/s${SUB}/${DIREXP}/anat_EPI_regTarget.nii.gz \
                    -float -final wsinc5
    done
done

## create intersession mask
for ((s=1;s<=SESS;s++)); do
    3dTstat -prefix ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/anat_EPI_tmp_motion_mean.nii.gz ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/func_motion1.nii.gz   
    3dAutomask -prefix ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/anat_EPI_mask_reg.nii.gz \
               -clfrac 0.5 -peels 3 -erode 1 \
               ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/anat_EPI_tmp_motion_mean.nii.gz
    rm ${DATA_DIR}/s${SUB}/${DIREXP}/SESS${s}/anat_EPI_tmp_motion_mean.nii.gz
done

3dcalc -prefix ${DATA_DIR}/s${SUB}/${DIREXP}/anat_EPI_mask.nii.gz \
       -a ${DATA_DIR}/s${SUB}/${DIREXP}/SESS1/anat_EPI_mask_reg.nii.gz \
       -b ${DATA_DIR}/s${SUB}/${DIREXP}/SESS2/anat_EPI_mask_reg.nii.gz \
       -expr 'and(a,b)'

done

## JDV Jul 12th 2013 ###########################################################
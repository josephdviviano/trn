## Dual Phase Encoded Stimuli Mean Generation w. AFNI & MATLAB
## Allows for multi-session registrations

## June 19th 2012, Joseph Viviano
## 
## SESS1 ... SESSn folders
##
# rename functional files in /SESSn

# This script:
# 1) Deobliques all PD scans and orients to RIA
# 2) Coregisters all PDs within session & Takes session mean
# 3) Coregisters all PD session means to session 1
# 4) Applies all transformations in one step using 3dAllineate at 2x resolution to deobliqued originals
# 5) Takes mean of all images at 2x resolution using an intersession mask

## SETTINGS ###################################################################
CORES=7       # number of CPU cores to use

fSES=1          # first session number
SESS=1          # number of sessions

RES=0.375      # upsampled ANAT resolution
FWHM=0.75      # smoothing full width half max (mm)
MEDIAN=3       # median filter num voxels for median (radius of cube)

INTERP=wsinc5

###############################################################################

## CREATE SESSION MEANS
## prepare each session for registration
for ((s=fSES;s<=SESS;s++)); do
    RUN=`ls -l SESS$s/raw/*PD* | wc -l`
    rm /tmp/list.txt; for ((i=1;i<=RUN;i++)); do echo $i >> /tmp/list.txt; done

    ## deoblique all PDs
    for ((i=1;i<=RUN;i++)); do
        3dWarp -deoblique -quintic -verb -prefix SESS$s/PD_ob$i.nii.gz SESS$s/raw/PD$i.nii.gz > SESS$s/mat_PD_tmp_deoblique_$i.1D
        sed '1 d' SESS$s/mat_PD_tmp_deoblique_$i.1D > SESS$s/mat_PD_tmp_deoblique_headerstrip_$i.1D
        cat_matvec -ONELINE SESS$s/mat_PD_tmp_deoblique_headerstrip_$i.1D > SESS$s/mat_PD_deoblique_$i.aff12.1D
        cat_matvec -ONELINE SESS$s/mat_PD_deoblique_$i.aff12.1D -I > SESS$s/mat_PD_deoblique_inv_$i.aff12.1D; done
    rm SESS$s/*mat_PD_tmp_deoblique*

    ## orient all files to RIA because BET is picky (RIA)
    cat /tmp/list.txt | parallel -P $CORES \
    3dAxialize -prefix SESS$s/PD_tmp_RIA{}.nii.gz -orient RIA SESS$s/PD_ob{}.nii.gz 

    ## skull strip each PD using stupid BET
    cat /tmp/list.txt | parallel -P $CORES \
    bet SESS$s/PD_tmp_RIA{}.nii.gz SESS$s/PD_tmp_mask_RIA{}.nii.gz -m -n -f 0.01 -g 0 -Z
    
    ## bring BET-derived mask back into the beautiful land of AFNI (RAI)
    cat /tmp/list.txt | parallel -P $CORES \
    3dAxialize -prefix SESS$s/PD_tmp_mask{}.nii.gz -orient RAI SESS$s/PD_tmp_mask_RIA{}_mask.nii.gz
    cat /tmp/list.txt | parallel -P $CORES \
    3dCalc -a SESS$s/PD_ob{}.nii.gz -b SESS$s/PD_tmp_mask{}.nii.gz -expr "'a*b'" -prefix SESS$s/PD_brain{}.nii.gz

    ## register each PD to the first scan of the session
    cd SESS$s
    for ((i=2;i<=RUN;i++)); do
    align_epi_anat.py -dset1 PD_brain1.nii.gz -dset2 PD_brain$i.nii.gz \
                      -dset2to1 -suffix PD_to_SESS \
                      -master_dset2 PD_brain1.nii.gz \
                      -dset1_strip None -dset2_strip None \
                      -volreg off -tshift off -deoblique off \
                      -big_move -partial_coronal \
                      -cost lpa

    rm PD_brain1.nii.gzPD_to_SESS_mat.aff12.1D
    mv PD_brain$i.nii.gzPD_to_SESS_mat.aff12.1D mat_PD_to_SESS_$i.aff12.1D
    mv PD_brain$i.nii.gzPD_to_SESS.nii.gz reg_PD_to_SESS_$i.nii.gz; done
    cp PD_brain1.nii.gz reg_PD_to_SESS_1.nii.gz
    cd ..

    ## create tmp session PDs and session masks
    3dMean -prefix anat_PD_tmp_mean_$s.nii.gz SESS$s/*reg_PD_to_SESS*
    3dAutomask -prefix SESS$s/anat_PD_tmp_mask_$s.nii.gz anat_PD_tmp_mean_$s.nii.gz

done

## REGISTER SESSION MEANS TO 1st SESSION
for ((s=2;s<=SESS;s++)); do
    align_epi_anat.py -dset1 anat_PD_tmp_mean_1.nii.gz \
                      -dset2 anat_PD_tmp_mean_$s.nii.gz \
                      -dset2to1 -suffix SESS_to_SESS \
                      -master_dset2 anat_PD_tmp_mean_1.nii.gz \
                      -dset1_strip None -dset2_strip None \
                      -volreg off -tshift off -deoblique off \
                      -giant_move -partial_coronal \
                      -cost lpa
    rm anat_PD_tmp_mean_1.nii.gzSESS_to_SESS_mat.aff12.1D
    mv anat_PD_tmp_mean_$s.nii.gzSESS_to_SESS_mat.aff12.1D mat_SESS_to_SESS_$s.aff12.1D
    mv anat_PD_tmp_mean_$s.nii.gzSESS_to_SESS.nii.gz reg_SESS_to_SESS_$s.nii.gz; done

## CREATE 2X RESOLUTION DUMMY FROM 1st SESSION
3dResample -dxyz $RES $RES $RES -rmode Li \
           -prefix anat_PD_resample_target.nii.gz \
           -inset  anat_PD_tmp_mean_1.nii.gz

## APPLY TRANSFORMATIONS TO ORIGINAL IMAGES at 2X RESOLUTION
# Apply to first run of first session
3dAllineate -prefix SESS1/proton_reg1.nii.gz \
            -input SESS1/raw/PD1.nii.gz \
            -1Dmatrix_apply SESS1/mat_PD_deoblique_inv_1.aff12.1D \
            -master anat_PD_resample_target.nii.gz \
            -float -final $INTERP
            
# Apply to rest of first session
RUN=`ls -l SESS1/raw/*PD* | wc -l`
mkdir SESS1/params1
for ((i=2;i<=RUN;i++)); do
    cat_matvec -ONELINE SESS1/mat_PD_deoblique_inv_$i.aff12.1D SESS1/mat_PD_to_SESS_$i.aff12.1D > SESS1/params1/reg$i.aff12.1D
    3dAllineate -prefix SESS1/proton_reg$i.nii.gz \
                -input SESS1/raw/PD$i.nii.gz \
                -1Dmatrix_apply SESS1/params1/reg$i.aff12.1D \
                -master anat_PD_resample_target.nii.gz \
                -float -final $INTERP
done

# Apply to first run of remaining sessions
for ((s=2;s<=SESS;s++)); do
    mkdir SESS$s/params$s
    cat_matvec -ONELINE SESS$s/mat_PD_deoblique_inv_1.aff12.1D mat_SESS_to_SESS_$s.aff12.1D > SESS$s/params$s/reg1.aff12.1D
    3dAllineate -prefix SESS$s/proton_reg1.nii.gz \
                -input SESS$s/raw/PD1.nii.gz \
                -1Dmatrix_apply SESS$s/params$s/reg1.aff12.1D \
                -master anat_PD_resample_target.nii.gz \
                -float -final $INTERP
done

# Apply to all remaining sessions
for ((s=2;s<=SESS;s++)); do
    RUN=`ls -l SESS$s/raw/*PD* | wc -l`
    for ((i=2;i<=RUN;i++)); do
    cat_matvec -ONELINE SESS$s/mat_PD_deoblique_inv_$i.aff12.1D SESS$s/mat_PD_to_SESS_$i.aff12.1D mat_SESS_to_SESS_$s.aff12.1D > SESS$s/params$s/reg$i.aff12.1D
    3dAllineate -prefix SESS$s/proton_reg$i.nii.gz \
                -input SESS$s/raw/PD$i.nii.gz \
                -1Dmatrix_apply SESS$s/params$s/reg$i.aff12.1D \
                -master anat_PD_resample_target.nii.gz \
                -float -final $INTERP
    done
done

## Smooth & Median Filtr all upsampled images
for ((s=fSES;s<=SESS;s++)); do
    RUN=`ls -l SESS$s/raw/*PD* | wc -l`
    rm /tmp/list.txt; for ((i=1;i<=RUN;i++)); do echo $i >> /tmp/list.txt; done
    
    cat /tmp/list.txt | parallel -P $CORES \
    3dBlurInMask -prefix SESS$s/proton_smooth{}.nii.gz \
                 -FWHM $FWHM \
                 -float \
                 -input SESS$s/proton_reg{}.nii.gz
    
    cat /tmp/list.txt | parallel -P $CORES \
    3dMedianFilter -prefix SESS$s/proton_median{}.nii.gz \
                   -irad $MEDIAN \
                   -iter 1 \
                   -verb \
                   SESS$s/proton_reg{}.nii.gz
    
#    cat /tmp/list.txt | parallel -P $CORES \
#    3dMedianFilter -prefix SESS$s/proton_median{}.nii.gz \
#                   -irad 5 \
#                   -iter 1 \
#                   -verb \
#                   SESS$s/proton_median_1stPass_{}.nii.gz
done

## Take Means (root mean square?)
3dMean -prefix anat_PD_smooth.nii.gz SESS1/*proton_smooth* # SESS2/*proton_reg* SESS3/*proton_reg* SESS4/*proton_reg*
3dMean -prefix anat_PD_median.nii.gz SESS1/*proton_median* # SESS2/*proton_reg* SESS3/*proton_reg* SESS4/*proton_reg*
3dMean -prefix anat_PD_reg.nii.gz SESS1/*proton_reg* 

## Deskull
3daxialize -prefix anat_PD_tmp_RIA.nii.gz -orient RIA anat_PD_reg.nii.gz
bet anat_PD_tmp_RIA.nii.gz PD_tmp_mask.nii.gz -m -n -f 0.25 -g -0.3 -Z
3daxialize -prefix PD_tmp_mask.nii.gz -orient RAI PD_tmp_mask_mask.nii.gz
3dcalc -prefix anat_PD_smooth_brain.nii.gz -a anat_PD_smooth.nii.gz -b PD_tmp_mask.nii.gz -expr 'a*b'
3dcalc -prefix anat_PD_median_brain.nii.gz -a anat_PD_median.nii.gz -b PD_tmp_mask.nii.gz -expr 'a*b'
3dcalc -prefix anat_PD_reg_brain.nii.gz -a anat_PD_reg.nii.gz -b PD_tmp_mask.nii.gz -expr 'a*b'
#3dAutobox -prefix anat_PD_smooth_brain.nii.gz -noclust -input anat_PD_tmp_smooth_brain.nii.gz # didnt work out to be the same number of voxels.
#3dAutobox -prefix anat_PD_median_brain.nii.gz -noclust -input anat_PD_tmp_median_brain.nii.gz # ...
rm *PD_tmp*

3daxialize -prefix anat_PD_tmp_RIA.nii.gz -orient RIA anat_PD_median.nii.gz
bet anat_PD_tmp_RIA.nii.gz PD_tmp_mask.nii.gz -m -n -f 0.25 -g -0.3 -Z
3daxialize -prefix PD_tmp_mask.nii.gz -orient RAI PD_tmp_mask_mask.nii.gz


## CREATE INTERSESSION MASK
#3dAutomask -prefix anat_PD_SESS_mask_1.nii.gz\
#                   anat_PD_tmp_mean_1.nii.gz
#
#for ((s=1;s<=SESS;s++));
#    do 3dAutomask -prefix anat_PD_SESS_mask_$s.nii.gz \
#                  reg_SESS_to_SESS_$s.nii.gz
#done
#
#3dcalc -prefix anat_PD_mask.nii.gz \
#       -a anat_PD_SESS_mask_1.nii.gz \
#       -b anat_PD_SESS_mask_2.nii.gz \
#       -c anat_PD_SESS_mask_3.nii.gz \
#       -d anat_PD_SESS_mask_4.nii.gz \
#       -expr 'and(a,b,c,d)'

#3dAllineate -prefix anat_PD_resample_mask.nii.gz \
#            -input anat_PD_mask.nii.gz \
#            -master anat_PD_resample_target.nii.gz \
#            -final NN -1Dparam_apply '1D: 12@0'\'

#3dCalc -prefix anat_PD_brain.nii.gz -a anat_PD.nii.gz -b anat_PD_resample_mask.nii.gz -expr 'a*b'

## The END

#3dAllineate -prefix anat_PD_SESS_mask_1.nii.gz \
#            -input SESS1/PD_tmp_mask1.nii.gz \
#            -master anat_PD_resample_target.nii.gz \
#            -final NN -1Dparam_apply '1D: 12@0'\'

#for ((s=2;s<=SESS;s++)); 
#    do 3dAllineate -prefix anat_PD_SESS_mask_$s.nii.gz \
#                   -input SESS$s/PD_tmp_mask$s.nii.gz \
#                   -1Dmatrix_apply mat_SESS_to_SESS_$s.aff12.1D \
#                   -master anat_PD_resample_target.nii.gz \
#                   -final NN; done



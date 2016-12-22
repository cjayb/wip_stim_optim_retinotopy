#!/bin/bash

SUBJECTS_DIR=/Users/cjb/projects/CFF_Retinotopy/fs_subjects_dir
subject=s18

surfreg --s $subject --t fsaverage_sym --lh
surfreg --s $subject --t fsaverage_sym --lh --xhemi

mri_surf2surf --srcsubject fsaverage_sym --srcsurfreg sphere.reg --trgsubject $subject --trgsurfreg fsaverage_sym.sphere.reg --hemi lh --sval angle-template-2.5.sym.mgh --tval ${SUBJECTS_DIR}/${subject}/label/lh.benson.angle.mgh

mri_surf2surf --srcsubject fsaverage_sym --srcsurfreg sphere.reg --trgsubject $subject --trgsurfreg fsaverage_sym.sphere.reg --hemi lh --sval eccen-template-2.5.sym.mgh --tval ${SUBJECTS_DIR}/${subject}/label/lh.benson.eccen.mgh

mri_surf2surf --srcsubject fsaverage_sym --srcsurfreg sphere.reg --trgsubject $subject --trgsurfreg fsaverage_sym.sphere.reg --hemi lh --sval areas-template-2.5.sym.mgh --tval ${SUBJECTS_DIR}/${subject}/label/lh.benson.areas.mgh

##
mri_surf2surf --srcsubject fsaverage_sym --srcsurfreg sphere.reg --trgsubject ${subject}/xhemi --trgsurfreg fsaverage_sym.sphere.reg --hemi lh --sval angle-template-2.5.sym.mgh --tval ${SUBJECTS_DIR}/${subject}/label/rh.benson.angle.mgh

mri_surf2surf --srcsubject fsaverage_sym --srcsurfreg sphere.reg --trgsubject ${subject}/xhemi --trgsurfreg fsaverage_sym.sphere.reg --hemi lh --sval eccen-template-2.5.sym.mgh --tval ${SUBJECTS_DIR}/${subject}/label/rh.benson.eccen.mgh

mri_surf2surf --srcsubject fsaverage_sym --srcsurfreg sphere.reg --trgsubject ${subject}/xhemi --trgsurfreg fsaverage_sym.sphere.reg --hemi lh --sval areas-template-2.5.sym.mgh --tval ${SUBJECTS_DIR}/${subject}/label/rh.benson.areas.mgh

#!/bin/bash

subject=030_WAH

mri_surf2surf --srcsubject fsaverage_sym --srcsurfreg sphere.reg --trgsubject $subject --trgsurfreg fsaverage_sym.sphere.reg --hemi lh --sval angle-template-2.5.sym.mgh --tval ${SUBJECTS_DIR}/${subject}/label/lh.benson.angle.mgh

mri_surf2surf --srcsubject fsaverage_sym --srcsurfreg sphere.reg --trgsubject $subject --trgsurfreg fsaverage_sym.sphere.reg --hemi lh --sval eccen-template-2.5.sym.mgh --tval ${SUBJECTS_DIR}/${subject}/label/lh.benson.eccen.mgh

mri_surf2surf --srcsubject fsaverage_sym --srcsurfreg sphere.reg --trgsubject $subject --trgsurfreg fsaverage_sym.sphere.reg --hemi lh --sval areas-template-2.5.sym.mgh --tval ${SUBJECTS_DIR}/${subject}/label/lh.benson.areas.mgh

## right hemisphere

# NBNBNB
# these were originally called with "--hemi lh"!
# Check that that's correct (should it not be "rh"?)!
mri_surf2surf --srcsubject fsaverage_sym --srcsurfreg sphere.reg --trgsubject ${subject}/xhemi --trgsurfreg fsaverage_sym.sphere.reg --hemi lh --sval angle-template-2.5.sym.mgh --tval ${SUBJECTS_DIR}/${subject}/label/rh.benson.angle.mgh

mri_surf2surf --srcsubject fsaverage_sym --srcsurfreg sphere.reg --trgsubject ${subject}/xhemi --trgsurfreg fsaverage_sym.sphere.reg --hemi lh --sval eccen-template-2.5.sym.mgh --tval ${SUBJECTS_DIR}/${subject}/label/rh.benson.eccen.mgh

mri_surf2surf --srcsubject fsaverage_sym --srcsurfreg sphere.reg --trgsubject ${subject}/xhemi --trgsurfreg fsaverage_sym.sphere.reg --hemi lh --sval areas-template-2.5.sym.mgh --tval ${SUBJECTS_DIR}/${subject}/label/rh.benson.areas.mgh

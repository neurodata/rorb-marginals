#!/bin/bash

### example neurodata.conf file here: <https://github.com/neurodata/ndex/blob/master/examples/neurodata.cfg.example>

### Separate call for the annotation layer

ndpull --config_file neurodata.conf --collection collman --experiment M247514_Rorb_1_Site3Align2_EM --channel m247514_Site3Annotation_MN_global --threads 4 --full_extent --outdir annotations/


chanArray=(synapsin_16bit PSD95_16bit)

export COLLECTION='collman'
export EXPERIMENT='M247514_Rorb_1_Site3Align2_LENS_Session1_CROP'
export BASE='data16'


for i in "${chanArray[@]}"
do 
  ndpull --config_file neurodata.conf --collection $COLLECTION --experiment $EXPERIMENT --channel $i --threads 4 --full_extent --outdir $BASE/$i/
done


### for 8-bit data

chanArray=(synapsin PSD95)

export COLLECTION='collman'
export EXPERIMENT='M247514_Rorb_1_Site3Align2_LENS_Session1_CROP'
export BASE='data'


for i in "${chanArray[@]}"
do 
  ndpull --config_file neurodata.conf --collection $COLLECTION --experiment $EXPERIMENT --channel $i --threads 4 --full_extent --outdir $BASE/$i/
done



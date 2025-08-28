#!/bin/bash

# Configuration
REMOTE_USER="mgorski"
LOCAL_PATH="/Users/mgorski/BADANIA/OCA/ocm_data"

# Separated pathr so that the -R (remote) option works well
REMOTE_PATH_ZB08="/work/apus/oca/fits/zb08/processed-ofp/targets"
REMOTE_PATH_WK06="/work/apus/oca/fits/wk06/processed-ofp/targets"

# Rsyncing matching files
rsync -avz -R \
      "$REMOTE_USER@araucaria.camk.edu.pl:$REMOTE_PATH_ZB08/./*/*/light-curve/*_diff_light_curve.txt" \
      "$LOCAL_PATH/zb08/"

rsync -avz -R \
      "$REMOTE_USER@araucaria.camk.edu.pl:$REMOTE_PATH_WK06/./*/*/light-curve/*_diff_light_curve.txt" \
      "$LOCAL_PATH/wk06/"
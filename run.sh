#!/bin/bash

ln -s ../data ./data 
matlab -nodisplay -nosoftwareopengl -r \
  "addpath(genpath('../code')); demo_2D;"
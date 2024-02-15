#!/bin/bash

echo "Z sample"
Rscript --vanilla alleleFreqModelZ_sample.R

echo "Z sim"
Rscript --vanilla alleleFreqModelZ_sim.R

echo "Z boot"
Rscript --vanilla alleleFreqModelZ_boot.R
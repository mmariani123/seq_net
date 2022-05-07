#!/usr/bin/env bash

#$ -q slipstream_queue@cn-t630
#$ -N run_rdnase
#$ -cwd
#$ -o /slipstream/home/mmariani/hhv6_ann/logs
#$ -j y

export PATH=/usr/bin/R:$PATH

Rscript \ 
/slipstream/home/mmariani/hhv6_ann/scripts/run_rdnase.R 

##$ -pe threads 16

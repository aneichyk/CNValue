#!/usr/bin/env Rscript

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2016 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Runs core TBRden function while generating no plots

#Load library
library(TBRden)

#Read args
args <- commandArgs(trailingOnly=T)

#Run TBRden
TBRden(args[1],args[2],OUTDIR=args[3],prefix=args[4],QQ=F,manhattan=F,gzip=F)
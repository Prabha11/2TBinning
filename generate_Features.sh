#!/bin/bash

name='simBG'
dir_path='sample_data/'
input_path=$dir_path
seqfile='sim.contig.nt'

# parameters
length_min=1000
length_max='1e99'

# OFDEG specific parameters
p_wordsize=4
p_stepsize=0.1
p_alpha=0.8
p_coverage=5

# -----------------------

perl generate_View1.pl $name $dir_path $input_path $seqfile $length_min $length_max $p_wordsize $p_stepsize $p_alpha $p_coverage
perl generate_View2.pl $name $dir_path $input_path $seqfile $length_min $length_max
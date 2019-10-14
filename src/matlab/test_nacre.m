% This is an example MATLAB script for quantifying nacre spots.
%
% Melissa R. McGuirl, Brown University. 2019.

% pfeffer example
input_dir = '../../data/sample_inputs/Out_nacre_default_1.mat';

% barcodes generated using:
%   python3 get_barcodes.py -i ../../data/sample_dist_mats/iriL_nacre_default_1_day56.txt -d 0 -o ../../data/sample_barcodes/iriL_nacre_default_1_day56
% and the input distance matrices were generated using processNacredata.m:
%   processNacredata('../../data/sample_inputs/', 'Out_nacre_default',  56, '../../data/sample_dist_mats', 1, '../../plots',  'nacre example')

PD_dir = '../../data/sample_barcodes/iriL_nacre_default_1_day56_dim0';

time_pt = 56;
pers_cutoff_iri = 100;
cell_type = 'I';

[num_spots, spot_size, roundness_score, center_stripe_rad, alignment_score,...
    xanS_mel_density, mel_xanS_density, mean_mel_space, var_mel_space, mean_melxanC_space, var_melxanC_space, ...
    melCV, mean_xanC_space, var_xanC_space]  = quantify_spots(input_dir, PD_dir, time_pt, pers_cutoff_iri, cell_type)

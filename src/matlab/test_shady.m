% This is an example MATLAB script for quantifying shady spots.
%
% Melissa R. McGuirl, Brown University. 2019.

% pfeffer example
input_dir = '../../data/sample_inputs/Out_shady_default_1.mat';

% barcodes generated using:
%   python3 get_barcodes.py -i ../../data/sample_dist_mats/melD_shady_default_1_day76.txt -d 0 -o ../../data/sample_barcodes/melD_shady_default_1_day76
% and the input distance matrices were generated using processShadyata.m:
%   processShadydata('../../data/sample_inputs/', 'Out_shady_default',  76, '../../data/sample_dist_mats', 1, '../../plots',  'shady example')

PD_dir = '../../data/sample_barcodes/melD_shady_default_1_day76_dim0';

time_pt = 76;
pers_cutoff_mel = 90;
cell_type = 'M';

[num_spots, spot_size, roundness_score, center_stripe_rad, alignment_score,...
    xanS_mel_density, mel_xanS_density, mean_mel_space, var_mel_space, mean_melxanC_space, var_melxanC_space, ...
    melCV, mean_xanC_space, var_xanC_space]  = quantify_spots(input_dir, PD_dir, time_pt, pers_cutoff_mel, cell_type)

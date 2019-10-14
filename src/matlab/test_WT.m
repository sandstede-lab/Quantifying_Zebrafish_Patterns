% WT example
input_dir = '../../data/sample_inputs/Out_WT_default_1.mat';

% barcodes generated using:
%   python3 get_barcodes.py -i ../../data/sample_dist_mats/melD_WT_default_1_day46.txt -d 1 -o ../../data/sample_barcodes//melD_WT_default_1_day46
%   python3 get_barcodes.py -i ../../data/sample_dist_mats/xanC_WT_default_1_day46.txt -d 1 -o ../../data/sample_barcodes/xanC_WT_default_1_day46
%   python3 get_barcodes.py -i ../../data/sample_dist_mats/xanS_WT_default_1_day46.txt -d 1 -o ../../data/sample_barcodes/xanS_WT_default_1_day46
% 
% and the input distance matrices were generated using processWTdata.m:
%   processWTdata('../../data/sample_inputs/', 'Out_WT_default', 46, '../../data/sample_dist_mats', 1, '../../plots',  'Wild-type example')



mel1_dir = '../../data/sample_barcodes/melD_WT_default_1_day46_dim1';
xanC1_dir =  '../../data/sample_barcodes/xanC_WT_default_1_day46_dim1';
xanS1_dir = '../../data/sample_barcodes/xanS_WT_default_1_day46_dim1';

time_pt = 46;

[num_stripes, num_Istripes, stripe_breaks, Istripe_breaks, dayOfNewStripes, max_stripe_separation,...
    max_Istripe_separation, avg_straightness, med_straightness,  ...
    mean_mel_space, var_mel_space, mean_xanC_space, var_xanC_space,...
    mean_xanS_space, var_xanS_space, mean_melxanC_space, var_melxanC_space, ...
    mean_melxanS_space, var_melxanS_space, melCV, xanS_mel_density, mel_xanS_density, ...
    iriLMel_density]  = quantify_stripes(input_dir, mel1_dir, xanC1_dir, xanS1_dir,  time_pt)



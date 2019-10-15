function  [num_stripes, num_Istripes, stripe_breaks, Istripe_breaks, dayOfNewStripes, max_stripe_separation, max_Istripe_separation, avg_straightness, med_straightness,  ...
     mean_mel_space, var_mel_space, mean_xanD_space, var_xanD_space,...
    mean_xanL_space, var_xanL_space, mean_melxanD_space, var_melxanD_space,  mean_melxanL_space, ...
    var_melxanL_space, melCV, xanL_mel_density, mel_xanL_density, iriLMel_density]  = quantify_stripes(cells_mel, cells_iriL, cells_xanD, cells_xanL,...
     mel1_dir, xanC1_dir, xanS1_dir, boundaryX, boundaryY, cellsXd_all, numXand_all, boundaryY_all)

% This is a MATLAB function whichs reads in files containing the
% dimension 0 and dimension 1 persistent homology data and outputs pattern
% statistics of zebrafish stripes.
%
% Inputs:cells_mel, cells_iriL, cells_xanD, and cells_xanL are coordinate data of cell locations at the time of interest; 
% mel1_dir, xanC1_dir, and xanS1_dir, are paths to the saved persistent homology data files (Ripser outputs); boundaryX
% is x-boundary of domain (right) at the time of interest; boundaryY is the y-boundary of the
% domain (top) at the time of interest; cellsXd_all is a matrix of coordinate data of cell
% locations of the X^d cells over time,  numXand_all is a vector containing
% the number of X^d cells over time, and boundary_all is a vector of
% boundaryY over time.
%
% Outputs: Complete collection of pattern statistcs for wild-type zebrafish
% based on topological data analysis, machine learning, and direct
% calculations.
%
% Author: Melissa R. McGuirl, Brown University. 2019.


% load in bar codes
bars_mel1 = importdata(mel1_dir);
bars_xanC1 = importdata(xanC1_dir);
bars_xanS1 = importdata(xanS1_dir);

%  get dimension 1 persistence; death-birth
mel1_pers = bars_mel1(:,2)-bars_mel1(:,1);
xanC1_pers = bars_xanC1(:,2)-bars_xanC1(:,1);
xanS1_pers = bars_xanS1(:,2)- bars_xanS1(:,1);

% define persistence thresholds based on cell-to-cell measurements
pers_cutoff = 200;

% compute dimension 1 betti numbers
b1_xanC = length(find(xanC1_pers > pers_cutoff  & bars_xanC1(:,1) < 80));
b1_xanS = length(find(xanS1_pers > pers_cutoff & bars_xanS1(:,1) < 100));
b1_mel = length(find(mel1_pers > pers_cutoff & bars_mel1(:,1) < 90));

num_stripes = b1_xanS; 
num_Istripes = b1_xanC;

if b1_xanS < 2 && b1_mel < 2
    stripe_breaks = 1;
else 
    stripe_breaks = 0;
end

if b1_xanC < 3 
    Istripe_breaks = 1;
else 
    Istripe_breaks = 0;
end

% use persistence of dim 1 features to get stripe widths, or maximum interstripe separation
% similarly for interstripe widths. exclude 'infinite' bar length by
% removing largest death point (first case). we see infinite bars in dim 1
% because of periodic boundary conditions.

xanC_widths = sort(xanC1_pers(xanC1_pers > pers_cutoff & bars_xanC1(:,1) < 80), 'descend');
xanC_widths = xanC_widths(2:end);
max_stripe_separation = 2*nanmedian(xanC_widths);

xanS_widths = sort(xanS1_pers(xanS1_pers > pers_cutoff  & bars_xanS1(:,1) < 100), 'descend');
xanS_widths = xanS_widths(2:end);
max_Istripe_separation = 2*nanmedian(xanS_widths);


% estimate first day when interstripes X1V and X1D begin to form
dayOfNewStripes = find_day_of_new_stripes(cellsXd_all, numXand_all, boundaryY_all);


% if M cells are present, compute local statistics
if ~isempty(cells_mel)
    
    % compute mean and variance of M-M nearest-neighbor distances
    [~, mean_mel_space, var_mel_space] = getPeriodicDistMats(cells_mel, boundaryX);
    
    % compute mean and variance of M-X^D nearest-neighbor distances
    D_MXD = pdist2(cells_mel, cells_xanD);
    indices = find(D_MXD > 120);
    D_MXD(indices) = nan;
    mean_melxanD_space = nanmean(min(D_MXD, [], 2));
    var_melxanD_space = nanvar(min(D_MXD, [], 2));
    
    % compute mean and variance of M-X^L nearest-neighbor distances
    D_MXL = pdist2(cells_mel, cells_xanD);
    indices = find(D_MXL > 120);
    D_MXL(indices) = nan;
    mean_melxanL_space = nanmean(min(D_MXL, [], 2));
    var_melxanL_space = nanvar(min(D_MXL, [], 2));
    
    % compute Mel CV.
    melCV =  100*(sqrt(var_mel_space)/mean_mel_space);
    
    % compute mean M-X^L density
    mel_xanL_density_temp = [];
    for k = 1:size(cells_mel,1)
        if  dist(cells_mel(k, 2), boundaryY) > 300 &&  dist(cells_mel(k, 1), boundaryX) > 300 && dist(cells_mel(k, 2), 0) > 300 &&  dist(cells_mel(k, 1), 0) > 300
            mel_xanL_density_temp(end+1) = length(find(dist(cells_mel(k,:),  cells_xanL') < 250));
        end
    end
    mel_xanL_density = min(quantile(mel_xanL_density_temp, 20));
    
else
    % If there are no M cells, set all local statistics to 0.
    mean_melxanL_space = 0;
    var_melxanL_space = 0;
    mean_melxanD_space = 0;
    var_melxanD_space = 0;
    melCV = 0;
    mean_mel_space = 0;
    var_mel_space = 0;
    mel_xanL_density = 0;
end

% Compute X^L-M densities. If there are no X^S cells, set density to 0.
if ~isempty(cells_xanL)
    xanL_mel_density_temp = [];
    for k = 1:size(cells_xanL,1)
        if  dist(cells_xanL(k, 2), boundaryY) > 300 &&  dist(cells_xanL(k, 1), boundaryX) > 300 && dist(cells_xanL(k, 2), 0) > 300 &&  dist(cells_xanL(k, 1), 0) > 300
            xanL_mel_density_temp(end+1) = length(find(dist(cells_xanL(k,:),   cells_mel') < 250));
        end
    end
    xanL_mel_density = min(quantile(xanL_mel_density_temp, 20));
else
    
    xanL_mel_density = 0;
end

% Compute I^L-M densities. If there are no I^L cells, set density to 0.
if ~isempty(cells_iriL)
    iri_mel_density = [];
    
    for k = 1:size(cells_iriL,1)
        if  dist(cells_iriL(k, 2), boundaryY) > 300 &&  dist(cells_iriL(k, 1), boundaryX) > 300 && dist(cells_iriL(k, 2), 0) > 300 &&  dist(cells_iriL(k, 1), 0) > 300
            iri_mel_density(end+1) = length(find(dist(cells_iriL(k,:),   cells_mel') < 200));
        end
    end
    iriLMel_density = min(quantile(iri_mel_density, 20));
else
    iriLMel_density = 0;
end


% compute mean and variance of X^L-X^L nearest-neighbor distances
if ~isempty(cells_xanL)
    [~, mean_xanL_space, var_xanL_space] = getPeriodicDistMats(cells_xanL, boundaryX);
else
    mean_xanL_space=0;
    var_xanL_space=0;
end


if ~isempty(cells_xanD)
    % compute mean and variance of X^C-X^C nearest-neighbor distances
    [~, mean_xanD_space, var_xanD_space] = getPeriodicDistMats(cells_xanD, boundaryX);
    
    % Use X^C cells to compute straightness measure of stripes 
    x_querys = 0:60:boundaryX;
    [top_cv, bottom_cv] = straightness_measure(cells_xanD, x_querys, b1_xanS, b1_mel, 'test', 0);
    
    avg_straightness = mean([top_cv; bottom_cv]);
    med_straightness = median([top_cv; bottom_cv]);
    
else
    % if there are no X^C cells, set statistics to 0.
    mean_xanD_space=0;
    var_xanD_space=0;
    avg_straightness = 0;
    med_straightness = 0;
end





end

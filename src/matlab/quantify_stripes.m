function  [num_stripes, num_Istripes, stripe_breaks, Istripe_breaks, dayOfNewStripes, max_stripe_separation, max_Istripe_separation, avg_straightness, med_straightness,  ...
     mean_mel_space, var_mel_space, mean_xanC_space, var_xanC_space,...
    mean_xanS_space, var_xanS_space, mean_melxanC_space, var_melxanC_space,  mean_melxanS_space, ...
    var_melxanS_space, melCV, xanS_mel_density, mel_xanS_density, iriLMel_density]  = quantify_stripes(output_dir, mel1_dir, xanC1_dir, xanS1_dir,  time_pt)

% This is a MATLAB function whichs reads in files containing the
% dimension 0 and dimension 1 persistent homology data and outputs pattern
% statistics of zebrafish stripes.
%
% Inputs: output_dir is the path to the original model output data; mel1_dir, xanC1_dir, 
% and xanS1_dir, are paths to the saved persistent homology data
% files (Ripser outputs); and time_point is the time/day at which
% patterns will be quantifies.
%
% Outputs: Complete collection of pattern statistcs for wild-type zebrafish
% based on topological data analysis, machine learning, and direct
% calculations.
%
% Author: Melissa R. McGuirl, Brown University. 2019.

% get cell positions
load(output_dir)

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
dayOfNewStripes = find_day_of_new_stripes(cellsXc, numXanc, boundaryY);

% extract cell coordinates by cell type, remove boundary cells
cutoff =  0.1*boundaryY(time_pt);
cells_mel = cellsM(find(cellsM(1:numMel(time_pt), 2,time_pt) > cutoff &  cellsM(1:numMel(time_pt), 2,time_pt) < boundaryY(time_pt) - cutoff), :, time_pt);
cells_xanC = cellsXc(find(cellsXc(1:numXanc(time_pt), 2,time_pt) > cutoff &  cellsXc(1:numXanc(time_pt), 2,time_pt) < boundaryY(time_pt) - cutoff),:, time_pt);
cells_xanS = cellsXsn(find(cellsXsn(1:numXansn(time_pt), 2,time_pt) > cutoff &  cellsXsn(1:numXansn(time_pt), 2,time_pt) < boundaryY(time_pt) - cutoff),:, time_pt);
cells_iriL = cellsIl(find(cellsIl(1:numIril(time_pt), 2,time_pt) > cutoff &  cellsIl(1:numIril(time_pt), 2,time_pt) < boundaryY(time_pt) - cutoff),:, time_pt);


% if M cells are present, compute local statistics
if ~isempty(cells_mel)
    
    % compute mean and variance of M-M nearest-neighbor distances
    [~, mean_mel_space, var_mel_space] = getPeriodicDistMats(cells_mel, boundaryX(time_pt));
    
    % compute mean and variance of M-X^C nearest-neighbor distances
    D_MXC = pdist2(cells_mel, cells_xanC);
    indices = find(D_MXC > 120);
    D_MXC(indices) = nan;
    mean_melxanC_space = nanmean(min(D_MXC, [], 2));
    var_melxanC_space = nanvar(min(D_MXC, [], 2));
    
    % compute mean and variance of M-X^S nearest-neighbor distances
    D_MXS = pdist2(cells_mel, cells_xanS);
    indices = find(D_MXS > 120);
    D_MXS(indices) = nan;
    mean_melxanS_space = nanmean(min(D_MXS, [], 2));
    var_melxanS_space = nanvar(min(D_MXS, [], 2));
    
    % compute Mel CV.
    melCV =  100*(sqrt(var_mel_space)/mean_mel_space);
    
    % compute mean M-xan^S density
    mel_xanS_density_temp = [];
    for k = 1:size(cells_mel,1)
        if  dist(cells_mel(k, 2), boundaryY(time_pt)) > 300 &&  dist(cells_mel(k, 1), boundaryX(time_pt)) > 300 && dist(cells_mel(k, 2), 0) > 300 &&  dist(cells_mel(k, 1), 0) > 300
            mel_xanS_density_temp(end+1) = length(find(dist(cells_mel(k,:),  cells_xanS') < 250));
        end
    end
    mel_xanS_density = min(quantile(mel_xanS_density_temp, 20));
    
else
    % If there are no M cells, set all local statistics to 0.
    mean_melxanS_space = 0;
    var_melxanS_space = 0;
    mean_melxanC_space = 0;
    var_melxanC_space = 0;
    melCV = 0;
    mean_mel_space = 0;
    var_mel_space = 0;
    mel_xanS_density = 0;
end

% Compute X^S-M densities. If there are no X^S cells, set density to 0.
if ~isempty(cells_xanS)
    xanS_mel_density_temp = [];
    for k = 1:size(cells_xanS,1)
        if  dist(cells_xanS(k, 2), boundaryY(time_pt)) > 300 &&  dist(cells_xanS(k, 1), boundaryX(time_pt)) > 300 && dist(cells_xanS(k, 2), 0) > 300 &&  dist(cells_xanS(k, 1), 0) > 300
            xanS_mel_density_temp(end+1) = length(find(dist(cells_xanS(k,:),   cells_mel') < 250));
        end
    end
    xanS_mel_density = min(quantile(xanS_mel_density_temp, 20));
else
    
    xanS_mel_density = 0;
end

% Compute I^L-M densities. If there are no I^L cells, set density to 0.
if ~isempty(cells_iriL)
    iri_mel_density = [];
    
    for k = 1:size(cells_iriL,1)
        if  dist(cells_iriL(k, 2), boundaryY(time_pt)) > 300 &&  dist(cells_iriL(k, 1), boundaryX(time_pt)) > 300 && dist(cells_iriL(k, 2), 0) > 300 &&  dist(cells_iriL(k, 1), 0) > 300
            iri_mel_density(end+1) = length(find(dist(cells_iriL(k,:),   cells_mel') < 200));
        end
    end
    iriLMel_density = min(quantile(iri_mel_density, 20));
else
    iriLMel_density = 0;
end


% compute mean and variance of X^S-X^S nearest-neighbor distances
if ~isempty(cells_xanS)
    [~, mean_xanS_space, var_xanS_space] = getPeriodicDistMats(cells_xanS, boundaryX(time_pt));
else
    mean_xanS_space=0;
    var_xanS_space=0;
end


if ~isempty(cells_xanC)
    % compute mean and variance of X^C-X^C nearest-neighbor distances
    [~, mean_xanC_space, var_xanC_space] = getPeriodicDistMats(cells_xanC, boundaryX(time_pt));
    
    % Use X^C cells to compute straightness measure of stripes 
    x_querys = 0:60:boundaryX(time_pt);
    [top_cv, bottom_cv] = straightness_measure(cells_xanC, x_querys, b1_xanS, b1_mel, 'test', 0);
    
    avg_straightness = mean([top_cv; bottom_cv]);
    med_straightness = median([top_cv; bottom_cv]);
    
else
    % if there are no X^C cells, set statistics to 0.
    mean_xanC_space=0;
    var_xanC_space=0;
    avg_straightness = 0;
    med_straightness = 0;
end





end

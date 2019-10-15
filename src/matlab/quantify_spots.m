function  [num_spots, spot_size, roundness_score, center_stripe_rad, alignment_score,...
    xanS_mel_density, mel_xanS_density, mean_mel_space, var_mel_space, mean_melxanC_space, var_melxanC_space, ...
    melCV, mean_xanC_space, var_xanC_space]  = quantify_spots(output_dir, PD_dir, time_pt, pers_cutoff, cell_type)

% This is a MATLAB function whichs reads in files containing the
% dimension 0 persistent homology data and outputs pattern
% statistics of a shady mutant.
%
% Inputs: output_dir is the path to the original model output data; OD_dir is a path
% to the saved persistent homology data files (Ripser outputs); time_point is the time/day at which
% patterns will be quantifies; pers_cutoff is the persistence cut for
% counting betti-0 numbers (based on cell-to-cell measurements); cell_type
% is the cell type used to compute spot statistics (options are 'M', 'I').
%
% Outputs: Complete collection of pattern statistcs for shady spots
% based on topological data analysis, machine learning, and direct
% calculations.
%
% Author: Melissa R. McGuirl, Brown University. 2019.

% load dimension 0 persistent homology data

% get cell positions
load(output_dir)

cutoff =  0.1*boundaryY(time_pt);
cells_mel = cellsM(find(cellsM(1:numMel(time_pt), 2,time_pt) > cutoff &  cellsM(1:numMel(time_pt), 2,time_pt) < boundaryY(time_pt) - cutoff), :, time_pt);
cells_iriL = cellsIl(find(cellsIl(1:numIril(time_pt), 2,time_pt) > cutoff &  cellsIl(1:numIril(time_pt), 2,time_pt) < boundaryY(time_pt) - cutoff),:, time_pt);
cells_xanC = cellsXc(find(cellsXc(1:numXanc(time_pt), 2,time_pt) > cutoff &  cellsXc(1:numXanc(time_pt), 2,time_pt) < boundaryY(time_pt) - cutoff),:, time_pt);
cells_xanS = cellsXsn(find(cellsXsn(1:numXansn(time_pt), 2,time_pt) > cutoff &  cellsXsn(1:numXansn(time_pt), 2,time_pt) < boundaryY(time_pt) - cutoff),:, time_pt);



bars_0 = importdata(PD_dir);

if isempty(bars_0)
    bars_0 = [0,0];
else
    bars_0(end,2) = bars_0(1,2);
    bars_0 = bars_0(2:end,:);
end

% dimension 0 so persistence = death since birth = 0
dim0_pers = bars_0(:,2);

% compute betii numbers
b0 = length(find(dim0_pers > pers_cutoff));
num_spots = b0; %number of spots = betti-0

if cell_type == 'M'
    cells = cells_mel;
elseif cell_type == 'I'
    cells = cells_iriL;
end

%cluster data to get spots, use b0 as number of clusters
if ~isempty(cells) && size(cells, 1) > 1 && b0 > 0
    if cell_type == 'M' % assume shady
        Z = linkage(cells, 'single',  '@my_dist_shady');
    else % assume pfeffer or nacre
        Z = linkage(cells, 'single',  '@my_dist_pfefnac');
    end
    c = cluster(Z,'maxclust',b0);
    cluster_size = zeros(b0,1);
    pca_ev_ratio = zeros(b0,1);
    centroids = zeros(b0,2);
    rad_spots = zeros(b0,1);
    
    for k = 1:b0
        cluster_size(k,1) = length(find(c == k)) ;
        if cluster_size(k,1) > 4 &&  boundaryX(time_pt) - max(cells(find(c==k),1)) > 500 && min(cells(find(c==k),1)) > 500 &&   boundaryY(time_pt) - max(cells(find(c==k),2)) > 500 && min(cells(find(c==k),2)) > 500
            [~, ~,latent] = pca(cells(find(c==k),:));
            pca_ev_ratio(k,1) = latent(1)/latent(2); % pca eigenvalue ratio
            centroids(k,:) = mean(cells(find(c == k),:)); % centroids of clusters
            rad_spots(k,1) = max(dist(cells(find(c== k),:), mean(cells(find(c== k),:))')); %radius of spots
        end
    end
else
    pca_ev_ratio = 0;
    cluster_size = 0;
end

spot_size = nanmedian(cluster_size); %spot size = median number of cells per cluster/spot
pca_ev_ratio = pca_ev_ratio( pca_ev_ratio ~= 0);
roundness_score = nanmedian(pca_ev_ratio);  %roundness score = median of pca eigenvalue ratios

% compute center stripe radius and alignment score
mid_stripe = boundaryY(time_pt)/2;
if ~isempty(cells)
    if size(cells, 1) > 1 && b0 > 4
        center_stripe_rad = min(dist(centroids, mid_stripe)) - median(rad_spots);
        D_centroids = squareform(pdist(centroids, 'chebychev'));
        nn = knnsearch(centroids, centroids, 'K', 2, 'Distance', 'chebychev');
        centroids_nnD = [];
        for k = 1:length(nn(:,1))
            centroids_nnD(k) = D_centroids(nn(k,1), nn(k,2));
        end
        alignment_score = std(centroids_nnD); %alignment score = standard deviation of NN distances (l-infinity)
        
    else
        % outlier case when there are not enough spots/cells
        center_stripe_rad = 0;
        alignment_score  = 1000000000000000000000;
    end
else
    center_stripe_rad = 0;
    alignment_score  = 1000000000000000000000;
    
end




if ~isempty(cells_xanS) && ~isempty(cells_mel)
    % compute X^S-M denstities
    xanS_mel_density_temp = [];
    for k = 1:size(cells_xanS,1)
        if  dist(cells_xanS(k, 2), boundaryY(time_pt)) > 300 &&  dist(cells_xanS(k, 1), boundaryX(time_pt)) > 300 && dist(cells_xanS(k, 2), 0) > 300 &&  dist(cells_xanS(k, 1), 0) > 300
            xanS_mel_density_temp(end+1) = length(find(dist(cells_xanS(k,:),   cells_mel') < 250));
        end
    end
    xanS_mel_density = min(quantile(xanS_mel_density_temp, 20));
    
    % compute M-X^S denstities
    mel_xanS_density_temp = [];
    for k = 1:size(cells_mel,1)
        if  dist(cells_mel(k, 2), boundaryY(time_pt)) > 300 &&  dist(cells_mel(k, 1), boundaryX(time_pt)) > 300 && dist(cells_mel(k, 2), 0) > 300 &&  dist(cells_mel(k, 1), 0) > 300
            mel_xanS_density_temp(end+1) = length(find(dist(cells_mel(k,:),  cells_xanS') < 250));
        end
    end
    mel_xanS_density = min(quantile(mel_xanS_density_temp, 20));
    
else
    xanS_mel_density = 0;
    mel_xanS_density = 0;
end

if ~isempty(cells_mel)
    
    % Compute M-M nearest neighbor distances
    [~, mean_mel_space, var_mel_space] = getPeriodicDistMats(cells_mel, boundaryX(time_pt));
    
    % Compute Mel CV
    melCV =  100*(sqrt(var_mel_space)/mean_mel_space);
    
else
    mean_mel_space = 0;
    var_mel_space = 0;
    melCV  = 0;
    
end

if  ~isempty(cells_mel) &&  ~isempty(cells_xanC)
    % Compute M-X^C nearest neighbor distances
    D_MXC = pdist2(cells_mel, cells_xanC);
    indices = find(D_MXC > 120); 
    D_MXC(indices) = nan;
    mean_melxanC_space = nanmean(min(D_MXC, [], 2));
    var_melxanC_space = nanvar(min(D_MXC, [], 2));
    
else
    mean_melxanC_space = 0;
    var_melxanC_space = 0;
end


% Compute X^C-X^C nearest neighbor distances
if ~isempty(cells_xanC)
    [~, mean_xanC_space, var_xanC_space] = getPeriodicDistMats(cells_xanC, boundaryX(time_pt));
else
    mean_xanC_space = 0;
    var_xanC_space = 0;
end


end

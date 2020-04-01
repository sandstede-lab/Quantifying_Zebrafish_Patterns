function  [num_spots, spot_size, roundness_score, center_stripe_rad, alignment_score,...
    xanL_mel_density, mel_xanL_density, mean_mel_space, var_mel_space, mean_melxanD_space, var_melxanD_space, ...
    melCV, mean_xanD_space, var_xanD_space]  = quantify_spots(cells_mel, cells_iriL, cells_xanD, cells_xanL, ...
    PD_dir, boundaryX, boundaryY, pers_cutoff, cell_type)

% This is a MATLAB function whichs reads in files containing the
% dimension 0 persistent homology data and outputs pattern
% statistics of a spotted mutant.
%
% Inputs: cells_mel, cells_iriL, cells_xanD, and cells_xanL are coordinate data of cell locations; PD_dir is a path
% to the saved persistent homology data files (Ripser outputs); boundaryX
% is x-boundary of domain (right); boundaryY is the y-boundary of the domain
% (top); pers_cutoff is the persistence cut for
% counting betti-0 numbers (based on cell-to-cell measurements); cell_type
% is the cell type used to compute spot statistics (options are 'M', 'I').
%
% Outputs: Complete collection of pattern statistcs for mutant spots
% based on topological data analysis, machine learning, and direct
% calculations.
%
% Author: Melissa R. McGuirl, Brown University. 2019.

% load dimension 0 persistent homology data
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
        if cluster_size(k,1) > 4 &&  boundaryX - max(cells(find(c==k),1)) > 500 && min(cells(find(c==k),1)) > 500 &&   boundaryY - max(cells(find(c==k),2)) > 500 && min(cells(find(c==k),2)) > 500
            [~, ~,latent] = pca(cells(find(c==k),:));
            pca_ev_ratio(k,1) = latent(1)/latent(2); % pca eigenvalue ratio for spot roundness
            centroids(k,:) = mean(cells(find(c == k),:)); % centroids of clusters
            rad_spots(k,1) = max(dist(cells(find(c== k),:), mean(cells(find(c== k),:))')); %radius of spots
        end
    end
else
    pca_ev_ratio = NaN;
    cluster_size = NaN;
end

spot_size = nanmedian(cluster_size); %spot size = median number of cells per cluster/spot
pca_ev_ratio = pca_ev_ratio( pca_ev_ratio ~= 0);
roundness_score = nanmedian(pca_ev_ratio);  %roundness score = median of pca eigenvalue ratios

% compute center stripe radius and alignment score
mid_stripe = boundaryY/2;
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
        center_stripe_rad = NaN;
        alignment_score  = NaN;
    end
else
    center_stripe_rad = NaN;
    alignment_score  = NaN;
    
end




if ~isempty(cells_xanL) && ~isempty(cells_mel)
    % compute X^L-M denstities
    xanL_mel_density_temp = [];
    for k = 1:size(cells_xanL,1)
        if  dist(cells_xanL(k, 2), boundaryY) > 300 &&  dist(cells_xanL(k, 1), boundaryX) > 300 && dist(cells_xanL(k, 2), 0) > 300 &&  dist(cells_xanL(k, 1), 0) > 300
            xanL_mel_density_temp(end+1) = length(find(dist(cells_xanL(k,:),   cells_mel') < 250));
        end
    end
    xanL_mel_density = min(quantile(xanL_mel_density_temp, 20));
    
    % compute M-X^L denstities
    mel_xanL_density_temp = [];
    for k = 1:size(cells_mel,1)
        if  dist(cells_mel(k, 2), boundaryY) > 300 &&  dist(cells_mel(k, 1), boundaryX) > 300 && dist(cells_mel(k, 2), 0) > 300 &&  dist(cells_mel(k, 1), 0) > 300
            mel_xanL_density_temp(end+1) = length(find(dist(cells_mel(k,:),  cells_xanL') < 250));
        end
    end
    mel_xanL_density = min(quantile(mel_xanL_density_temp, 20));
    
else
    xanL_mel_density = NaN;
    mel_xanL_density = NaN;
end

if ~isempty(cells_mel)
    
    % Compute M-M nearest neighbor distances
    [~, mean_mel_space, var_mel_space] = getPeriodicDistMats(cells_mel, boundaryX);
    
    % Compute Mel CV
    melCV =  100*(sqrt(var_mel_space)/mean_mel_space);
    
else
    mean_mel_space = NaN;
    var_mel_space = NaN;
    melCV  = NaN;
    
end

if  ~isempty(cells_mel) &&  ~isempty(cells_xanD)
    % Compute M-X^D nearest neighbor distances
    D_MXC = pdist2(cells_mel, cells_xanD);
    indices = find(D_MXC > 120); 
    D_MXC(indices) = nan;
    mean_melxanD_space = nanmean(min(D_MXC, [], 2));
    var_melxanD_space = nanvar(min(D_MXC, [], 2));
    
else
    mean_melxanD_space = NaN;
    var_melxanD_space = NaN;
end


% Compute X^D-X^D nearest neighbor distances
if ~isempty(cells_xanD)
    [~, mean_xanD_space, var_xanD_space] = getPeriodicDistMats(cells_xanD, boundaryX);
else
    mean_xanD_space = NaN;
    var_xanD_space = NaN;
end


end

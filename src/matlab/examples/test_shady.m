% This is an example MATLAB script for quantifying shady spots.
%
% Melissa R. McGuirl, Brown University. 2019.

% shady  example

clear all; close all;

addpath ../util
addpath ../

% input MAT file from model simulation 
input_dir = '../../../data/sample_inputs/Out_shady_default_1.mat';

% load in data to get cell positions
load(input_dir)
time_pt = 76;

% extract cell coordinates by cell type, remove boundary cells
cutoff =  0.1*boundaryY(time_pt);
cells_mel = cellsM(find(cellsM(1:numMel(time_pt), 2,time_pt) > cutoff &  cellsM(1:numMel(time_pt), 2,time_pt) < boundaryY(time_pt) - cutoff), :, time_pt);
cells_iriL = cellsIl(find(cellsIl(1:numIril(time_pt), 2,time_pt) > cutoff &  cellsIl(1:numIril(time_pt), 2,time_pt) < boundaryY(time_pt) - cutoff),:, time_pt);
cells_iriD = cellsId(find(cellsId(1:numIrid(time_pt), 2,time_pt) > cutoff &  cellsId(1:numIrid(time_pt), 2,time_pt) < boundaryY(time_pt) - cutoff),:, time_pt);
cells_xanD = cellsXc(find(cellsXc(1:numXanc(time_pt), 2,time_pt) > cutoff &  cellsXc(1:numXanc(time_pt), 2,time_pt) < boundaryY(time_pt) - cutoff),:, time_pt);
cells_xanL = cellsXsn(find(cellsXsn(1:numXansn(time_pt), 2,time_pt) > cutoff &  cellsXsn(1:numXansn(time_pt), 2,time_pt) < boundaryY(time_pt) - cutoff),:, time_pt);

%get distances
D_mel= getPeriodicDistMats(cells_mel, boundaryX(time_pt));
     
% Save data to text files 
savename = 'Shady_Test';
path2outputs = '../../../data/sample_dist_mats';
dlmwrite([path2outputs '/melD_' savename '_day' num2str(time_pt) '.txt'], D_mel, 'precision',3)
    
    
plotOn = 1; %set to 0 to shuf off plotting 
% plot if specified
if plotOn == 1
    f = figure(1);
    rectangle('Position',[-50,-50,boundaryX(time_pt)+100,boundaryY(time_pt)+100],'FaceColor',[255/255 239/255 213/255],'EdgeColor',[255/255 239/255 213/255],...
        'LineWidth',1)
    hold on;
    plot(cells_mel(:,1), cells_mel(:,2), 'o', 'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize', 3)
    plot(cells_xanD(:,1), cells_xanD(:,2), '*', 'MarkerEdgeColor',[255/255  185/255 15/255],'MarkerFaceColor',[255/255  185/255 15/255],'MarkerSize', 3)
    plot(cells_xanL(:,1), cells_xanL(:,2), 'yx', 'MarkerEdgeColor',[255/255  185/255 15/255],'MarkerSize', 2)
    plot(cells_iriL(:,1), cells_iriL(:,2), 'v','MarkerEdgeColor',[58/255 95/255 205/255],'MarkerSize', 4)
    plot(cells_iriD(:,1), cells_iriD(:,2), 's','MarkerEdgeColor',[192/255 192/255 192/255],'MarkerSize', 4)

    title('Shady example')
    xlim([-250 5000])
    set(gca,'XTick',[])
    set(gca,'YTick',[])
end
    
% define boundary 
boundaryX = boundaryX(time_pt);
boundaryY = boundaryY(time_pt);

% instructions to get persistent homology data 
fprintf('Next use Python to run Ripser to generate barcodes via: \n')
fprintf('cd ../../python \n')
fprintf('python3 get_barcodes.py -i ../../data/sample_dist_mats/melD_Shady_Test_day76.txt -d 1 -o ../../data/sample_barcodes/melD_Shady_Test_day76 \n')
pause; 


PD_dir = '../../../data/sample_barcodes/melD_shady_Test_day76_dim0';


%persistence threshold and cell-type for spot quantification
pers_cutoff_mel = 90;
cell_type = 'M';

% quantify patterns 
[num_spots, spot_size, roundness_score, center_stripe_rad, alignment_score,...
    xanL_mel_density, mel_xanL_density, mean_mel_space, var_mel_space, mean_melxanD_space, var_melxanD_space, ...
    melCV, mean_xanD_space, var_xanD_space]  = quantify_spots(cells_mel, cells_iriL, cells_xanD, cells_xanL, ...
    PD_dir, boundaryX, boundaryY, pers_cutoff_mel, cell_type)

% This is an example MATLAB script for quantifying zebrafish stripes.
%
% Melissa R. McGuirl, Brown University. 2019.

% WT example

clear all; close all;

addpath ../util
addpath ../


% input MAT file from model simulation 
input_dir = '../../../data/sample_inputs/Out_WT_default_1.mat';

% load in data to get cell positions
load(input_dir)
time_pt = 46;

% extract cell coordinates by cell type, remove boundary cells
cutoff =  0.1*boundaryY(time_pt);
cells_mel = cellsM(find(cellsM(1:numMel(time_pt), 2,time_pt) > cutoff &  cellsM(1:numMel(time_pt), 2,time_pt) < boundaryY(time_pt) - cutoff), :, time_pt);
cells_xanD = cellsXc(find(cellsXc(1:numXanc(time_pt), 2,time_pt) > cutoff &  cellsXc(1:numXanc(time_pt), 2,time_pt) < boundaryY(time_pt) - cutoff),:, time_pt);
cells_xanL = cellsXsn(find(cellsXsn(1:numXansn(time_pt), 2,time_pt) > cutoff &  cellsXsn(1:numXansn(time_pt), 2,time_pt) < boundaryY(time_pt) - cutoff),:, time_pt);
cells_iriL = cellsIl(find(cellsIl(1:numIril(time_pt), 2,time_pt) > cutoff &  cellsIl(1:numIril(time_pt), 2,time_pt) < boundaryY(time_pt) - cutoff),:, time_pt);
cells_iriD = cellsId(find(cellsId(1:numIrid(time_pt), 2,time_pt) > cutoff &  cellsId(1:numIrid(time_pt), 2,time_pt) < boundaryY(time_pt) - cutoff),:, time_pt);
 
%get distance matrices 
[D_mel, ~, ~] = getPeriodicDistMats(cells_mel, boundaryX(time_pt));
[D_xanD, ~, ~] = getPeriodicDistMats(cells_xanD, boundaryX(time_pt));
[D_xanL, ~, ~] = getPeriodicDistMats(cells_xanL, boundaryX(time_pt));

 % Save data to text files 
 savename = 'WT_Test';
 path2outputs = '../../../data/sample_dist_mats';
 dlmwrite([path2outputs '/melD_' savename '_day' num2str(time_pt) '.txt'], D_mel, 'precision',3)
 dlmwrite([path2outputs '/xanD_' savename '_day' num2str(time_pt) '.txt'], D_xanD, 'precision',3)
 dlmwrite([path2outputs '/xanL_' savename '_day' num2str(time_pt) '.txt'], D_xanL, 'precision',3)
 
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

    title('Wild type example')
    xlim([-250 5000])
    set(gca,'XTick',[])
    set(gca,'YTick',[])
end
    
% boundary of domain
boundaryY_all = boundaryY; 
boundaryY = boundaryY(time_pt);
boundaryX = boundaryX(time_pt);

% time series of X^d cells for identifying pattern events (when 2nd and 3rd
% interstripes form)
cellsXd_all = cellsXc;
numXand_all = numXanc;

% instructions to get persistent homology data 
fprintf('Next use Python to run Ripser to generate barcodes via: \n')
fprintf('cd ../../python \n')
fprintf('python3 get_barcodes.py -i ../../data/sample_dist_mats/melD_WT_Test_day46.txt -d 1 -o ../../data/sample_barcodes/melD_WT_Test_day46 \n')
fprintf('python3 get_barcodes.py -i ../../data/sample_dist_mats/xanD_WT_Test_day46.txt -d 1 -o ../../data/sample_barcodes/xanD_WT_Test_day46 \n')
fprintf('python3 get_barcodes.py -i ../../data/sample_dist_mats/xanL_WT_Test_day46.txt -d 1 -o ../../data/sample_barcodes/xanL_WT_Test_day46 \n')
pause; 

mel1_dir = '../../../data/sample_barcodes/melD_WT_Test_day46_dim1';
xanD1_dir =  '../../../data/sample_barcodes/xanD_WT_Test_day46_dim1';
xanL1_dir = '../../../data/sample_barcodes/xanL_WT_Test_day46_dim1';

% quantify patterns 
[num_stripes, num_Istripes, stripe_breaks, Istripe_breaks, dayOfNewStripes, stripe_width, Istripe_width, avg_straightness, med_straightness,  ...
     mean_mel_space, var_mel_space, mean_xanD_space, var_xanD_space,...
    mean_xanL_space, var_xanL_space, mean_melxanD_space, var_melxanD_space,  mean_melxanL_space, ...
    var_melxanL_space, melCV, xanL_mel_density, mel_xanL_density, iriLMel_density]  = quantify_stripes(cells_mel, cells_iriL, cells_xanD, cells_xanL,...
     mel1_dir, xanD1_dir, xanL1_dir, boundaryX, boundaryY, cellsXd_all, numXand_all, boundaryY_all)



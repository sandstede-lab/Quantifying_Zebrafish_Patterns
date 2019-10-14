function [] = processNacredata(path2inputs, name_scheme, time_pt, path2outputs, plotOn, plot_path, plot_title)

% This is a MATLAB function whichs reads in input .mat files containing the
% black and yellow cell coordinates over a time frame, and extracts
% the data from the specified time point. One distance matrix with periodic
% boundary conditions is saved as a text file and the spots are plotted and saved if plot = 1.
%
% Inputs: path2inputs is path to input data folder; name_scheme specifies
% how inputs of interest are labeled; time_pt is the time/day at which
% patterns will be quantifies; path2outputs is path to output folder;
% plotOn is a flag for plotting patterns (= 1 to plot, =0 to not plot);
% plot_path is path to save plots and plot_title is a string variable for
% the title of your plot. 
%
% Outputs: No variable outputs. This function saves distance matrices as
% text files and saves plots if specified. 
% 
% This function depends on getPeriodicDistMats.m.
%
% Author: Melissa R. McGuirl, Brown University. 2019.

DIR = dir([path2inputs '/' name_scheme '*']);

for i = 1:length(DIR)
    

    % load data
    load([path2inputs '/' DIR(i).name])
    
    % extract cells time point data, remove boundaries as 10% 
    cutoff =  0.1*boundaryY(time_pt);
    
    cells_xanC = cellsXc(find(cellsXc(1:numXanc(time_pt), 2,time_pt) > cutoff &  cellsXc(1:numXanc(time_pt), 2,time_pt) < boundaryY(time_pt) - cutoff),:, time_pt);
    cells_xanS = cellsXsn(find(cellsXsn(1:numXansn(time_pt), 2,time_pt) > cutoff &  cellsXsn(1:numXansn(time_pt), 2,time_pt) < boundaryY(time_pt) - cutoff),:, time_pt);
    cells_iriD = cellsId(find(cellsId(1:numIrid(time_pt), 2,time_pt) > cutoff &  cellsId(1:numIrid(time_pt), 2,time_pt) < boundaryY(time_pt) - cutoff),:, time_pt);
    cells_iriL = cellsIl(find(cellsIl(1:numIril(time_pt), 2,time_pt) > cutoff &  cellsIl(1:numIril(time_pt), 2,time_pt) < boundaryY(time_pt) - cutoff),:, time_pt);
    
    
    %get distances
    D_iriL = getPeriodicDistMats(cells_iriL, boundaryX(time_pt));
    % Save data
    dlmwrite([path2outputs '/iriL_' DIR(i).name(5:end-4) '_day' num2str(time_pt) '.txt'], D_iriL, 'precision',3)

    % plot if specified
    if plotOn == 1
        %plot black and yellow cells
        f = figure(1);
        rectangle('Position',[-50,-50,boundaryX(time_pt)+100,boundaryY(time_pt)+100],'FaceColor',[255/255 239/255 213/255],'EdgeColor',[255/255 239/255 213/255],...
            'LineWidth',1)
        hold on;
        plot(cells_xanC(:,1), cells_xanC(:,2), '*', 'MarkerEdgeColor',[255/255  185/255 15/255],'MarkerFaceColor',[255/255  185/255 15/255],'MarkerSize', 3)
        plot(cells_xanS(:,1), cells_xanS(:,2), 'yx', 'MarkerEdgeColor',[255/255  185/255 15/255],'MarkerSize', 2)
        plot(cells_iriL(:,1), cells_iriL(:,2), 'v','MarkerEdgeColor',[58/255 95/255 205/255],'MarkerSize', 4)
        plot(cells_iriD(:,1), cells_iriD(:,2), 's','MarkerEdgeColor',[192/255 192/255 192/255],'MarkerSize', 4)
         
        title(plot_title)
        xlim([-250 5000])
        set(gca,'XTick',[])
        set(gca,'YTick',[])
        saveas(f, strcat(plot_path, '/', DIR(i).name(5:end-4), '.jpg'))
        close all
    end
    
    
end


end






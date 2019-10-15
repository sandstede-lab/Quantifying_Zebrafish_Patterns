function [top_cv, bottom_cv] = straightness_measure(cells_xanD, x_querys, xanL_betti_1, mel_betti_1, plotName, plotOn)

% This is a MATLAB function that uses topological data analysis and
% clustering to compute a straightness measure of stripe patterns.
%
% Inputs: cells_xanD are the coordinates of X^d cells in the pattern of
% interest; x_querys are the x-coordinates of the grid points used for the
% numerical computation of the arc length distance of the stripe boundary; xanL_betti_1
% is the first betti number of X^l cells; mel_betti_1 is the first betti
% number of the M cells; plotName is save name of plot; and plotOn is a flag for plotting patterns (= 1 to plot, =0 to not plot)
%
% Outputs: straightness measures for the top and bottom boundary of each stripe 
%
% Author: Melissa R. McGuirl, Brown University. 2019.


% Apply single linkage clusters to X^d cells with number of clusters n,
% set to number of expected stripes (3) minus the number of stripe breaks,
% as indicated by the betti numbers.

Z_xanD = linkage(cells_xanD, 'single', '@my_dist_WT');

if xanL_betti_1 >= 2 || mel_betti_1 >= 2
    n = 3;
elseif xanL_betti_1 == 1 || mel_betti_1 == 1
    n = 2;
else
    n = 1;
end
c_xanD = cluster(Z_xanD,'maxclust', n);

% create variables to store straightness measure statistics
top_cv = zeros(n,1);
bottom_cv = zeros(n,1);

% loop through each stripe
for j = 1:n
    
    % find cells which correspond to stripe cluster
    my_cluster = cells_xanD(find(c_xanD == j),:);
    
    % get x-coordinates of boundary cells
    top_bd_x = zeros(length(x_querys)-1, 1);
    bottom_bd_x = zeros(length(x_querys)-1, 1);
    
    % create variable for y-coordinates of boundary cells
    top_bd_y = zeros(length(x_querys)-1, 1);
    bottom_bd_y = zeros(length(x_querys)-1, 1);
    
    % For each grid point, find the maximum and minimum y-coordinates over
    % all cells in the grid. These correspond to the boundary cells.
    
    flagOn = [];
    startflagOn = 0;
    for k = 1:length(x_querys)-1
        
        if ~isempty(find(my_cluster(:,1) >= x_querys(k) &  my_cluster(:,1) < x_querys(k + 1)))
            top_bd_x(k,1) = max(my_cluster(find(my_cluster(:,1) >= x_querys(k) &  my_cluster(:,1) < x_querys(k + 1)), 1));
            top_bd_y(k,1) = max(my_cluster(find(my_cluster(:,1) >= x_querys(k) &  my_cluster(:,1) < x_querys(k + 1)), 2));
            bottom_bd_x(k,1) = min(my_cluster(find(my_cluster(:,1) >= x_querys(k) &  my_cluster(:,1) < x_querys(k + 1)), 1));
            bottom_bd_y(k,1) = min(my_cluster(find(my_cluster(:,1) >= x_querys(k) &  my_cluster(:,1) < x_querys(k + 1)), 2));
            
            % flag to handle breaks
            if k > 1
                if abs(top_bd_y(k,1) - top_bd_y(k-1,1)) > 300 || abs(bottom_bd_y(k,1) - bottom_bd_y(k-1,1)) > 300
                    flagOn(end + 1) = k;
                end
            end
            
        else % if there are no cells in the current grid point, take a straight line distance between adjacent grids
            if k > 1
                top_bd_x(k,1) = top_bd_x(k-1,1) + min(diff(x_querys));
                top_bd_y(k,1) = top_bd_y(k-1,1);
                bottom_bd_x(k,1) = bottom_bd_x(k-1,1) + min(diff(x_querys));
                bottom_bd_y(k,1) = bottom_bd_y(k-1,1);
                flagOn(end + 1) = k;
            else
                startflagOn = k;
            end
            
        end
    end
    
    % handles cases for stripe breaks
    if startflagOn > 0
        top_bd_x = top_bd_x(startflagOn +1:end);
        bottom_bd_x = bottom_bd_x(startflagOn +1:end);
        top_bd_y = top_bd_y(startflagOn +1:end);
        bottom_bd_y = bottom_bd_y(startflagOn +1:end);
    end
    
    if isempty(flagOn)
        % is no breaks, compute straighness measure directly
        top_cv(j,1) = (sum(sqrt(diff(top_bd_x).^2 + diff(top_bd_y).^2))/max(x_querys) - 1)*100;
        bottom_cv(j,1)=(sum(sqrt(diff(bottom_bd_x).^2 + diff(bottom_bd_y).^2))/max(x_querys) - 1)*100;
        
        if plotOn == 1
            % plot  stripes with measure
            f = figure(1);
            hold on
            scatter(my_cluster(:,1), my_cluster(:,2))
            plot(top_bd_x, top_bd_y, '-k', 'lineWidth', 2)
            plot(bottom_bd_x, bottom_bd_y,  '-k','lineWidth', 2)
            text(50 + x_querys(end-1),top_bd_y(end) + 5, [num2str(top_cv(j,1)) '%'], 'FontSize',24)
            text(50 + x_querys(end-1),bottom_bd_y(end) - 5, [num2str(bottom_cv(j,1)) '%'], 'FontSize',24)
        end
        
    else
        
        %if a stripe break occurred, compute arc length distances and stripe straightness measure
        %separately for each component of broken stripe.
        
        top_bd_x_1 = top_bd_x(1:min(flagOn)-1);
        top_bd_y_1 = top_bd_y(1:min(flagOn)-1);
        bottom_bd_x_1 =  bottom_bd_x(1:min(flagOn)-1);
        bottom_bd_y_1 = bottom_bd_y(1:min(flagOn)-1);
        
        top_bd_x_2 = top_bd_x(max(flagOn)+1:end);
        top_bd_y_2 = top_bd_y(max(flagOn)+1:end);
        bottom_bd_x_2 =  bottom_bd_x(max(flagOn)+1:end);
        bottom_bd_y_2 = bottom_bd_y(max(flagOn)+1:end);
        
        top_cv(j,1) = ((sum(sqrt(diff(top_bd_x_1).^2 + diff(top_bd_y_1).^2)) + ...
            sum(sqrt(diff(top_bd_x_2).^2 + diff(top_bd_y_2).^2)) +...
            abs(x_querys(max(flagOn)) - x_querys(min(flagOn))))/max(x_querys) - 1)*100;
        
        bottom_cv(j,1) = ((sum(sqrt(diff(bottom_bd_x_1).^2 + diff(bottom_bd_y_1).^2)) + ...
            sum(sqrt(diff(bottom_bd_x_2).^2 + diff(bottom_bd_y_2).^2)) +...
            abs(x_querys(max(flagOn)) - x_querys(min(flagOn))))/max(x_querys) - 1)*100;
        
        
        if plotOn == 1
            f = figure(1);
            % plot  stripes with measure
            hold on
            scatter(my_cluster(:,1), my_cluster(:,2))
            plot(top_bd_x_1, top_bd_y_1, '-k', 'lineWidth', 2)
            plot(bottom_bd_x_1, bottom_bd_y_1,  '-k','lineWidth', 2)
            plot(top_bd_x_2, top_bd_y_2, '-k', 'lineWidth', 2)
            plot(bottom_bd_x_2, bottom_bd_y_2,  '-k','lineWidth', 2)
            text(50 + x_querys(end-1),top_bd_y(end) + 5, [num2str(top_cv(j,1)) '%'], 'FontSize',24)
            text(50 + x_querys(end-1),bottom_bd_y(end) - 5, [num2str(bottom_cv(j,1)) '%'], 'FontSize',24)
        end
        
        
        
    end
    
    
    
end


% Save plot
if plotOn == 1
    f = figure(1);
    set(gca,'fontsize',24);
    xlim([0 max(x_querys) + 2000])
    saveas(f, plotName)
    clf(f)
end






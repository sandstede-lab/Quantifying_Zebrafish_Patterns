function [D, mean_dists, var_dists]  = getPeriodicDistMats(coords, boundaryX)
% This is a MATLAB function which takes as input a collection of x-y coordinates and
% computes the pairwise distances between all coordinates, with periodic boundary conditions 
% imposed in the x direction. The periodicity condition is impossed with
% the boundaryX input parameter, which specifies the uppper limit of the
% x-coordinates. We assume the lower limit of the x-coordinates is 0. 
%
% Inputs: coords is list of x-y coordinates; boundaryX is the right
% x-direction boundary of the domain.
%
% Outputs: D is a distance matrix between input coordinates, mean_dists is
% the mean nearest-neighbor distance, and var_dists is the variance of the
% nearest-neighbor distances.
%
% Author: Melissa R. McGuirl, Brown University. 2019.

    mydistfun = @(x,y)(sqrt((x(:,2)-y(:,2)).^2 + min((x(:,1)-y(:,1)).^2, (boundaryX -abs(x(:,1)-y(:,1))).^2)));
    D = pdist2(coords,coords,mydistfun);
    
    %compute mean and variance of distance to closest cell
    mean_dists = mean(min(eye(size(D,1),size(D,2))*1000 + D, [], 2));
    var_dists = var(min(eye(size(D,1),size(D,2))*1000 + D, [], 2));
    
    % Ensures output is nonempty 
    if isempty(D)
        D = 0;
    end 
 
end


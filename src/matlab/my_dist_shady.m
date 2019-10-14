function D = my_dist_shady(x, y)
%boundaryX = 4850; %3710 for WT, 4850 for shady, 4090 for nacre and pfef
boundaryX = 4850;
D = (sqrt((x(:,2)-y(:,2)).^2 + min((x(:,1)-y(:,1)).^2, ( boundaryX -abs(x(:,1)-y(:,1))).^2)));
end 
function dayOfNewStripes = find_day_of_new_stripes(cellsXc, numXanc, boundaryY) 

% This is a MATLAB function that estimates the first day when interstripes X1V and X1D begin to form
%
% Inputs: cellsXc are the coordinates of X^C cells in the pattern of
% interest; numXanc is a vector of the number of xan^C cells at each time
% point, boundaryY is vector of the y-boundaries at each time point.
%
% Outputs: dayOfNewStripes is the estimated day when interstripes X1V and X1D begin to form
%
% Author: Melissa R. McGuirl, Brown University. 2019.


% assume interstripes X1V and X1D begin to form somewhere between 32 and 62
% dpf. Start at day 31 as IC. 
times = [10:1:41];
lowerbounds = zeros(length(times), 1);
upperbounds = zeros(length(times), 1);
dayOfNewStripes = 0;
count = 1;
% continue until we see new interstripes forming, or until we exceed 62 dpf 
while dayOfNewStripes == 0 && count <= length(times)
    time_current = times(count);
    cutoff =  0.1*boundaryY( time_current);
    cells_xanC = cellsXc(find(cellsXc(1:numXanc( time_current), 2, time_current) > cutoff &  cellsXc(1:numXanc( time_current), 2, time_current) < boundaryY( time_current) - cutoff),:,  time_current);
    if ~isempty(cells_xanC)
        lowerbounds(count,1) = min(cells_xanC(:,2));
        upperbounds(count,1) = max(cells_xanC(:,2));
        if count > 1 
            % If new cells begins to form > 200 mum away from center
            % stripe, assume that new interstripes are starting to form.
            if lowerbounds(count, 1) - lowerbounds(count-1 , 1) > 200 || upperbounds(count,1) - upperbounds(count-1, 1) > 200
                dayOfNewStripes  =  time_current;
            end
        end
    end
    count = count + 1;
end



end 
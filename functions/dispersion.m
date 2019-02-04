function disp = dispersion(x,y)
% This works out the dispersion (average distance to the center of a bunch
% of coordinates). x and y are the list of x and y coordinates for each
% point.

% Work out the center of the points
meanx = mean(x);
meany = mean(y);

% calculate distances to center
for i = 1:length(x);
    distance(i) = sqrt( (x(i) - meanx)^2 + (y(i) - meany)^2 );
end

disp = mean(distance);
end
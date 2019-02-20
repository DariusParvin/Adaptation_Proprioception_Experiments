function disp = dispersion2(P1x,P1y,P2x,P2y)
% This works out the dispersion (average distance to the center of a bunch
% of coordinates). x and y are the list of x and y coordinates for each
% point.
% The '2' refers to using two different centroids for block 1 and block 2


% calculate distances to center for trials part 1
for i = 1:length(P1x);
    distances1(i) = sqrt( (P1x(i) - mean(P1x))^2 + (P1y(i) - mean(P1y))^2 );
end

% calculate distances to center for trials part 2
for i = 1:length(P2x);
    distances2(i) = sqrt( (P2x(i) - mean(P2x))^2 + (P2y(i) - mean(P2y))^2 );
end

disp = mean([distances1 distances2]);
end
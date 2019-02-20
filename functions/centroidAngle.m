function cAngle = centroidAngle(x,y)
% This works out the angle between the hand position, start position (-160mm relative to hand), and the centroid of points

P0 = [0,-160]; % Start position 
P1 = [0,0]; % True hand position
P2 = [mean(x),mean(y)]; % centroid of proprioceptive estimates

% Calculate the angle between P0-P1 and P0-P2
cAngle = atan2(abs(det([P2-P0;P1-P0])),dot(P2-P0,P1-P0))*180/pi;
cAngle = cAngle * -sign(mean(x)); % Positive is counterclockwise 
end
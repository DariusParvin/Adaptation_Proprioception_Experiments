function sem_y = sem(x)

% if min(size(x)) > 1
%     DP: error - 'please only put in a one dimensional array, I dont know if this works for 2d array'
% end 

std_x = nanstd(x); 
n_x = sum(~isnan(x));

sem_y = std_x./sqrt(n_x);

end
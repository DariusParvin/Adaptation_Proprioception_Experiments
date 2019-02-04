function plot_correlation(data, xVar, yVar)

figure; hold on

% Get x data
xVar = {xVar};
xData = data.(xVar{:});

% Get y data
yVar = {yVar};
yData = data.(yVar{:});

plot(xData, yData ,'.','markersize',20)

[rho,pval] = corrcoef(xData,yData);
r = rho(1,2);
p = pval(1,2);

xlabel(xVar)
ylabel(yVar)

str = sprintf('%s vs %s, r = %.2f, p = %.2f',xVar{:}, yVar{:}, r, p);
title(str)

end
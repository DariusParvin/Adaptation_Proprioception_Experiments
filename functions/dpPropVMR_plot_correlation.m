function dpPropVMR_plot_correlation(data, xVar, yVar)
hold on

% Get x and y data
xVar = {xVar};
xData = data.(xVar{:});
yVar = {yVar};
yData = data.(yVar{:});

plot(xData, yData ,'.','markersize',20)

% show r and p value
[rho,pval] = corrcoef(xData,yData);
str = sprintf('r = %.3f, p = %.3f',rho(1,2),pval(1,2) );

if pval(1,2) < 0.05  % Bold text if significant
    text(0.07,0.93,str,'Units','normalized','FontSize',8,'FontWeight','Bold')
else
    text(0.07,0.93,str,'Units','normalized','FontSize',8)
end

xlabel(xVar)
ylabel(yVar)
str = sprintf('%s vs %s,',xVar{:}, yVar{:});
title(str)
set(gca,'FontSize',8)
end
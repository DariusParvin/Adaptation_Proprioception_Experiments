function dpPlotBar(x,dataPoints)

% Plots a single bar on x with datapoints y.
% Plots SEM bars and individual points

% x = x position of bar
% y = y data points

y = nanmean(dataPoints);
err = sem(dataPoints);
bar(x, y, .75); hold on; % plot bar
plot([x x], [y-err y+err],'linewidth',2,'color','k');   % error bar

jitter = randn(length(dataPoints), 1)./50;
offset = 0.1;
plot(repmat(x,length(y),1) + jitter + offset  , dataPoints,...
    'ok','markersize',4,'markeredgecolor','k','markerfacecolor','w') % plot points

end
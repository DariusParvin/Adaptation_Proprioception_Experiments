function dpPlotBar2(x, dataPoints, faceColor, plot_individuals)
% Plots a single bar on x with datapoints y and SEM bars

y = nanmean(dataPoints);
err = sem(dataPoints);
bar(x, y, .75, 'FaceColor',faceColor); hold on; % plot bar
plot([x x], [y-err y+err],'linewidth',2,'color','k');   % error bar


if nargin > 3
    jitter = randn(length(dataPoints), 1)./50;
    offset = 0.15;
    plot(repmat(x,length(y),1) + jitter + offset  , dataPoints,...
        'ok','markersize',3,'markeredgecolor','k','markerfacecolor','w') % plot points
end

end
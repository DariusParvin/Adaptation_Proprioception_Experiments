function dpPlotBar2(x,dataPoints,faceColor)
% Plots a single bar on x with datapoints y and SEM bars

y = nanmean(dataPoints);
err = sem(dataPoints);
bar(x, y, .75, 'FaceColor',faceColor); hold on; % plot bar
plot([x x], [y-err y+err],'linewidth',2,'color','k');   % error bar

end
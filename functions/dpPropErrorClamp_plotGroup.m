function dpPropVMR_plotGroup(data, varToPlot, blocksToPlot, axisLim, lineCol)

% data = data in table format
% varToPlot = y variable
% blocksToPlot = x variable - reaching blocks ( 'RB' ) or proprioceptive blocks ( 'PB' )
% axisLim = axis limits in the format [xmin xmax ymin ymax]
% lineCol = line color e.g. 'r'


% Get x data
blocksToPlot = {blocksToPlot};
blockData = data.(blocksToPlot{:}) ;

% Get y data
varToPlot = {varToPlot};
varData = data.(varToPlot{:});

% Plot hand angle for reaching blocks
for bi = 1:max( blockData )
    block_trials=[]; y=[]; sem_y=[];
    idx = ( blockData == bi ); %  block index
    block_trials = unique(data.TN(idx)); %  block trials
    
    for i = 1:length(block_trials); % loop through blocks
        trial_idx = data.TN==block_trials(i);
        y(i) = nanmean( varData(trial_idx ) );
        sem_y(i) = sem( varData(trial_idx ) );
    end
    
    % Plot hand angles
    shadedErrorBar(block_trials', y, sem_y, lineCol)
           
end


% Figure reference lines
axis(axisLim);

% Draw block breaks
drawline1(0, 'dir', 'horz', 'linestyle', '-'); 
blockLines = [9.5 36.5 72.5 108.5 144.5 180.5 270.5 360.5 396.5 486.5 522.5 612.5 648.5 738.5];
drawline1(blockLines , 'dir', 'vert', 'linestyle', ':');  % Draws line for blocks

% Shade the no feedback trials
no_fb_base =patch([9.5 36.5 36.5 9.5],[min(ylim) min(ylim) max(ylim) max(ylim)],zeros(1,4));
set(no_fb_base,'facecolor',[0 0 0]); set(no_fb_base,'edgealpha',0);
alpha(no_fb_base,0.1)

end



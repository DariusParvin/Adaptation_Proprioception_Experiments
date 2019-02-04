%% Set directory
clearvars
addpath(genpath('..'))
% save figures here
figDir = '../Figures/';

%% Load and clean Data
load('pre-processed_data/AlissaPropEC_trials.mat')

K_subj_names = {'PropVMR_602A__S53A', 'PropVMR_602A__S54A', ... % begin day 1
    'PropVMR_602A__S55A', 'PropVMR_602A__S56A', ...
    'PropVMR_602A__S58A', ...
    'PropVMR_602A__S59A', 'PropVMR_602A__S61A', ...
    'PropVMR_602A__S62A', 'PropVMR_602A__S63A', ...
    'PropVMR_602A__S64A', 'PropVMR_602A__S65A', ...
    'PropVMR_602A__S66A', 'PropVMR_602A__S67A', ...
    'PropVMR_602A__S69A', 'PropVMR_602A__S70A', ...
    'PropVMR_602C__S36A', 'PropVMR_602C__S46A', ...
    'PropVMR_602C__S37A', 'PropVMR_602C__S38A', ...
    'PropVMR_602C__S39A', 'PropVMR_602C__S40A', ...
    'PropVMR_602C__S41A', 'PropVMR_602C__S42A', ...
    'PropVMR_602C__S43A', 'PropVMR_602C__S44A',...
    'PropVMR_602C__S45A', 'PropVMR_602C__S47A', ...
    'PropVMR_602C__S48A', 'PropVMR_602C__S50A', ...
    'PropVMR_602C__S51A', ... % end of day 1
    'PropVMR_602A__S53B', 'PropVMR_602A__S54B', ... % begin day 2
    'PropVMR_602A__S55B', 'PropVMR_602A__S56B', ...
    'PropVMR_602A__S58B', ...
    'PropVMR_602A__S59B', 'PropVMR_602A__S61B', ...
    'PropVMR_602A__S62B', 'PropVMR_602A__S63B', ...
    'PropVMR_602A__S64B', 'PropVMR_602A__S65B', ...
    'PropVMR_602A__S66B', 'PropVMR_602A__S67B', ...
    'PropVMR_602A__S69B', 'PropVMR_602A__S70B', ...
    'PropVMR_602C__S36B', 'PropVMR_602C__S46B', ...
    'PropVMR_602C__S37B', 'PropVMR_602C__S38B', ...
    'PropVMR_602C__S39B', 'PropVMR_602C__S40B', ...
    'PropVMR_602C__S41B', 'PropVMR_602C__S42B', ...
    'PropVMR_602C__S43B', 'PropVMR_602C__S44B',...
    'PropVMR_602C__S45B', 'PropVMR_602C__S47B', ...
    'PropVMR_602C__S48B', 'PropVMR_602C__S50B', ...
    'PropVMR_602C__S51B'
    };

% Check that # of subject names == # of subjects in datafile
if length(K_subj_names)/2 ~= length(unique(T.SN))
    error_please_copy_in_subject_names_T1
end

prop_vars = {'FC_bias_X', 'FC_bias_Y', 'prop_theta'};

T.hand = T.hand_theta;
remove_vars = {'hand_theta','hand_theta_maxv','hand_theta_maxradv','handMaxRadExt','hand_theta_50'};
T(:, T.Properties.VariableNames(remove_vars)) = [];

% T1 REMOVE OUTLIERS
T1 = T;
outlier_idx = abs(T1.hand) > 90 ; % Remove trials greater than x degrees
fprintf('Outlier trials removed: %d \n' , sum(outlier_idx))
T1.hand(outlier_idx, 1) = nan; % Flip trials .*(-1)


% T2 FLIP CCW SUBJECTS TO CW
T2 = T1;
flip_idx = T2.rot_cond > 0; % CW condition index
T2.hand(flip_idx, 1) = T1.hand(flip_idx, 1).*(-1); % Flip trials .*(-1)

% flip proprioceptive related variables
T2.prop_theta(flip_idx) = T2.prop_theta(flip_idx).*(-1);
T2.FC_X(flip_idx) = T2.FC_X(flip_idx).*(-1);
T2.HL_X(flip_idx) = T2.HL_X(flip_idx).*(-1);
T2.FC_bias_X(flip_idx) = T2.FC_bias_X(flip_idx).*(-1);
T2.PropLocX(flip_idx) = T2.PropLocX(flip_idx).*(-1);
T2.ti(flip_idx) = T2.ti(flip_idx).*(-1) + 180;
T2.PropTestAng(flip_idx) = T2.PropTestAng(flip_idx).*(-1) + 180;


% T3 BASELINE SUBTRACTION
T3 = T2;
baseCN = 3:4; %%%% Baseline cycles to subtract %%%% CHECK THIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

base_idx = T3.CN >= min(baseCN) & T3.CN <= max(baseCN); % index of baseline cycles
base_mean = varfun(@nanmean,T2(base_idx ,:),'GroupingVariables',{'SN','day_num','ti'},'OutputFormat','table');



for SN = unique(T3.SN)'
    for di = unique(T3.day_num)'
        for ti = unique(T3.ti(~isnan(T3.ti)))' % subtract baseline for each target
            trial_idx = (T3.SN==SN & T.day_num == di & T3.ti==ti);
            base_idx = (base_mean.SN==SN & base_mean.day_num==di & base_mean.ti==ti);
            T3.hand(trial_idx) = T2.hand(trial_idx) - base_mean.nanmean_hand(base_idx);
        end
    end
end

prop_baseCN = 11:14; %%%% Proprioceptive Baseline cycles to subtract
% prop_baseCN = 5:8; %%%% Proprioceptive Baseline cycles to subtract
prop_base_idx = T3.CN >= min(prop_baseCN) & T3.CN <= max(prop_baseCN); % index of baseline cycles
prop_base_mean = varfun(@nanmean,T2(prop_base_idx ,:),'GroupingVariables',{'SN','day_num','PropTestAng'},'OutputFormat','table');

for SN = unique(T3.SN)'
    for di = unique(T3.day_num)'
        for prop_ti = unique(T3.PropTestAng(~isnan(T3.PropTestAng)))' % subtract baseline for each target
            for vi = 1:length(prop_vars) % loop over hand angle columns
                trial_idx = (T3.SN==SN & T3.day_num == di & T3.PropTestAng==prop_ti);
                prop_base_idx = (prop_base_mean.SN==SN & prop_base_mean.day_num == di &  prop_base_mean.PropTestAng==prop_ti);
                T3.(prop_vars{vi})(trial_idx) = T2.(prop_vars{vi})(trial_idx) - prop_base_mean.(strcat('nanmean_',prop_vars{vi}))(prop_base_idx);
            end
        end
    end
end

% Figure lines for block divisions
K_lines_all_trials = [5.5 10.5 20.5 70.5 160.5 165.5 195.5 235.5 265.5 305.5 310.5 340.5 380.5 385.5 415.5 420.5 470.5];

%% Individual subject, every trial
close all;
E = T3;
E.PB(isnan(E.PB)) = 0; % Cheap hack so that the split function works
day = 2; % day 1 or 2

subjs = unique(E.SN);

% for i = 1:length(subjs)   % for all subj
for i = 1 % subj to plot
    indx = (E.SN == subjs(i) & E.day_num == day);

    figure('units','centimeters','pos',[1 5 20 20]);hold on;
    subplot(2,2,1:2); hold on;
    str = sprintf('Hand angle and Proprioceptive estimates for subj %d, day %d', subjs(i), day);
    title(str);
    % Hand theta
    x1 = E.TN(indx);
    y1 = E.hand(indx);
    scatter(x1,y1,10,'filled');

    % Proprioceptive estimate (as an angle relative to tgt)
    x2 = E.TN(indx);
    y2 = E.prop_theta(indx);
    scatter(x2,y2,10,'filled');
    % Reference lines
    drawline1(0, 'dir', 'horz', 'linestyle', '-'); % Draws line at target
    drawline1(K_lines_all_trials, 'dir', 'vert', 'linestyle', ':');  % Draws line for blocks
    % Shade the no feedback trials
    no_fb_base =patch([0 10.5 10.5 0],[min(ylim) min(ylim) max(ylim) max(ylim)],zeros(1,4));
    set(no_fb_base,'facecolor',[0 0 0]); set(no_fb_base,'edgealpha',0);
    alpha(no_fb_base,0.1)
    % Fig labels
    xlabel('Trial'); ylabel('Hand Angle/Proprioceptive estimate (º)')

    subplot(2,2,3); title('Proprioceptive estimates'); hold on;
    % All proprioceptive estimates in absolute space
    x3 = E.FC_bias_X(indx);
    y3 = E.FC_bias_Y(indx);
    %     x3 = E.FC_X(indx);
    %     y3 = E.FC_Y(indx);
    Prop_block = E.PB(indx);
    scatterplot(x3,y3,'split',Prop_block,'leg','auto' );
    % Fig labels
    axis('square');
    xlabel('mm'); ylabel('mm')


    subplot(2,2,4); title('Mean Proprioceptive Bias'); hold on;
    % Mean Proprioceptive estimate per block
    x3=[]; y3=[];
    PBs = [1:6]; % Proprioceptive block number
    for bi = 1:length(PBs)
        x3(bi,1) = nanmean( E.FC_bias_X(indx & E.PB == PBs(bi)) );
        y3(bi,1) = nanmean( E.FC_bias_Y(indx & E.PB == PBs(bi)) );
    end
    cmap = jet(length(x3)); % Color map
    scatter(x3,y3,40,cmap,'filled');
    text(x3,y3,{'   bnf','   bf','   c1','   c2','   c3','   c4'}) % Label points

    % Fig labels
    scatter(0,0,80,'k.') % Target reference (0,0)
    %     axis([-30 30 -30 30]);
    axis('square');
    xlabel('mm'); ylabel('mm')
end

% print(sprintf('%sE1egSubj_%s',figDir,date),'-painters','-dpdf')
print(sprintf('%sE1egSubj_%s',figDir,date),'-painters','-djpeg')

%% Group hand angle and proprioceptive bias
day = 1;
E = T3(T3.day_num  == day, :); %day 1 or 2
figure; hold on;
set(gcf,'units','centimeters','pos',[5 5 20 10]);

dpPropVMR_plotGroup(E, 'hand', 'RB', [0 max(E.TN) -15 35], 'b')
dpPropVMR_plotGroup(E, 'prop_theta', 'PB', [0 max(E.TN) -15 35], 'r')

xlabel('Trial'); ylabel('Hand Angle/Proprioceptive estimate (º)')

% print(sprintf('%E1_Group_Hand_%s',figDir,date),'-painters','-dpdf')
print(sprintf('%sE1_Group_day%d_%s',figDir,day,date),'-painters','-djpeg')

%% Subject reliability table
clearvars -except T* K* figDir*
E = T3;
subj = unique(E.SN);
days = unique(E.day_num);

sdi = 0; % subject and day index key
for si = 1:length(subj)
    for di = 1:length(days)
        
        subj_day_idx = (E.SN==subj(si) & E.day_num==days(di));
        sdi = sdi+1;
        
        SN(sdi,1) = subj(si);
        Day(sdi,1) = days(di);
        
        rotCond(sdi,1) = mean(E.rot_cond(subj_day_idx,1));
        
        % Adaptation dependent variables
        afterEffect(sdi,1) = nanmean( E.hand( subj_day_idx & E.CN >= 77  & E.CN <= 77 ) ); % last aftereffect
%         afterEffect(sdi,1) = nanmean( E.hand( subj_day_idx & (E.CN == 62 | E.CN == 77 | E.CN == 84 ) ) ); % last aftereffect

        % Proprioceptive shift
        shift_idx = subj_day_idx & E.PB > 2 ;
        base_idx = subj_day_idx & E.PB == 1 ;
        propShiftX(sdi,1) = nanmean(E.FC_bias_X(shift_idx)) - nanmean(E.FC_bias_X(base_idx));
        
        propShiftTheta(sdi,1) = nanmean(E.prop_theta(shift_idx)) - nanmean(E.prop_theta(base_idx));
        
        % Proprioceptive dispersion
        disp_idx = subj_day_idx & E.PB==1;
        dispBlock1(sdi,1) = dispersion( E.FC_bias_X(disp_idx), E.FC_bias_Y(disp_idx));
        
        disp_idx=[];
        disp_idx = subj_day_idx & ~isnan(E.PB);
        dispAll(sdi,1) = dispersion( E.FC_bias_X(disp_idx), E.FC_bias_Y(disp_idx));
    end
end
% Summary table
reliabilityMatrix = table(afterEffect, propShiftX, propShiftTheta, dispBlock1, dispAll);

day1Matrix = reliabilityMatrix(Day==1,:);
day2Matrix = reliabilityMatrix(Day==2,:);

%% Subject summary table
clearvars -except T* K* summary* day* figDir* 
E = T3;
subj = unique(E.SN);

for si = 1:length(subj)
    
    SN(si,1) = subj(si);
    rotCond(si,1) = mean(E.rot_cond(E.SN==subj(si),1));
    
    % Adaptation dependent variables
    earlyLearning(si,1) = nanmean( E.hand( E.SN==subj(si) & E.CN >= (21 + 1)  & E.CN <= (21 + 4) ) ); % clamp cycles 2 - 4 (smaller clamp cycles than usual because of generalization)
    afterEffect(si,1) = nanmean( E.hand( E.SN==subj(si) & E.CN >= 77  & E.CN <= 77 ) ); % last aftereffect
    
    % Proprioceptive shift
    shift_idx = E.SN==subj(si) & E.PB > 2 ;
    base_idx = E.SN==subj(si) & E.PB == 1 ;
    propShiftX(si,1) = nanmean(E.FC_bias_X(shift_idx)) - nanmean(E.FC_bias_X(base_idx));
    
    propShiftTheta(si,1) = nanmean(E.prop_theta(shift_idx)) - nanmean(E.prop_theta(base_idx));
    
    % Proprioceptive dispersion
    disp_idx = E.SN==subj(si) & E.PB==1;
    dispBlock1(si,1) = dispersion( E.FC_bias_X(disp_idx), E.FC_bias_Y(disp_idx));
    
    disp_idx=[];
    disp_idx = E.SN==subj(si) & ~isnan(E.PB);
    dispAll(si,1) = dispersion( E.FC_bias_X(disp_idx), E.FC_bias_Y(disp_idx));
end
% Summary table
summaryMatrix = table(earlyLearning, afterEffect, propShiftX, propShiftTheta, dispBlock1, dispAll);

subjCond = table(SN, rotCond);
summaryBar = [summaryMatrix subjCond];

%% RELIABILITY CORRELATION MATRIX
clearvars -except T* K* summary* day* figDir* 
figure; set(gcf,'units','centimeters','pos',[5 5 30 20]);

varnames = day1Matrix.Properties.VariableNames;
for vi = 1:length(varnames) % loop over hand angle columns
    subplot(2, ceil(length(varnames)/2), vi);
    
    x = day1Matrix.(varnames{vi});
    y = day2Matrix.(varnames{vi});
    AX = plot(x,y,'.','markersize',8);
    
    % Axes and reference line
    axis('square')
    min_ax = min([xlim ylim]);
    max_ax = max([xlim ylim]);
    axis( [min_ax*0.9 max_ax*1.1 min_ax*0.9 max_ax*1.1 ]);    
    rline = refline(1,0);    
    rline.Color = [.5 .5 .5];
    rline.LineStyle = '-';
    
    [rho,pval] = corrcoef(x,y);
    str = sprintf('r = %.3f, p = %.3f',rho(1,2),pval(1,2) );
    %     text(xlim(1),ylim(2),str)
    text(0.05,0.9,str,'Units','normalized','FontSize',10)
    title(sprintf('%s',varnames{vi}) );
    xlabel('Day 1');
    ylabel('Day 2');
    set(gca, 'FontName', 'Arial', 'FontSize',12);
end

% print(sprintf('%E1_reliability_%s',figDir,date),'-painters','-dpdf')
print(sprintf('%sE1_reliability_%s',figDir,date),'-painters','-djpeg')
%% CORRELATION MATRIX
figure; set(gcf,'units','centimeters','pos',[5 5 20 20]);
[S,AX,BigAx,H,HAx] = plotmatrix(summaryMatrix{:,:});
[rho,pval] = corrcoef(summaryMatrix{:,:});

almostSigPlots = find(pval<0.1);
sigPlots = find(pval<0.05);

varnames = summaryMatrix.Properties.VariableNames;
for vi = 1:length(varnames) % loop over hand angle columns
    AX(vi,1).YLabel.String = varnames{vi};
    AX(length(varnames),vi).XLabel.String = varnames{vi};
end

for sigi = 1:length(almostSigPlots)
    S(almostSigPlots(sigi)).Color = [1 0.5 0];
end

for sigi = 1:length(sigPlots)
    S(sigPlots(sigi)).Color = 'r';
end

% print(sprintf('%E1_corr_Matrix_%s',figDir,date),'-painters','-dpdf')
print(sprintf('%sE1_corr_Matrix_%s',figDir,date),'-painters','-djpeg')

%% Plot specific correlation
dpPropVMR_plot_correlation(summaryMatrix, 'dispAll', 'afterEffect')

set(gcf,'units','centimeters','pos',[5 5 15 15]);

    % Axes and reference line
    axis('square')
    min_ax = min([xlim ylim]);
    max_ax = max([xlim ylim]);
    axis( [min_ax*0.9 max_ax*1.1 min_ax*0.9 max_ax*1.1 ]);    
    rline = refline(1,0);    
    rline.Color = [.5 .5 .5];
    rline.LineStyle = '-';

% print(sprintf('%E1_disp_vs_asymp_%s',figDir,date),'-painters','-dpdf')
print(sprintf('%sE1_disp_vs_asymp_%s',figDir,date),'-painters','-djpeg')

%% Plot specific correlation
dpPropVMR_plot_correlation(summaryMatrix, 'dispAll', 'propShiftTheta')

set(gcf,'units','centimeters','pos',[5 5 15 15]);

    % Axes and reference line
    axis('square')
    min_ax = min([xlim ylim]);
    max_ax = max([xlim ylim]);
    axis( [min_ax*0.9 max_ax*1.1 min_ax*0.9 max_ax*1.1 ]);    
    rline = refline(1,0);    
    rline.Color = [.5 .5 .5];
    rline.LineStyle = '-';

% print(sprintf('%sE1_disp_vs_shift_%s',figDir,date),'-painters','-dpdf')
% print(sprintf('%sE1_disp_vs_shift_%s',figDir,date),'-painters','-djpeg')

%% BAR GRAPHS Split CCW and CW
figure
% Plot bars

barVarNames = summaryBar.Properties.VariableNames;
xpos = [1 2 4 5 7 8 10 11 13 14 16 17];
for vi = 1:6
    dataPoints = summaryBar.(barVarNames{vi})(summaryBar.rotCond < 0);
    dpPlotBar(xpos(vi*2), dataPoints );
    
    dataPoints = summaryBar.(barVarNames{vi})(summaryBar.rotCond > 0);
    dpPlotBar(xpos(vi*2 - 1), dataPoints );
end

% Aesthetics
xticks = [1.5 4.5 7.5 10.5 13.5 16.5];
% set(gca,'xTick',xticks,'YTick',yticks,'xticklabel',{'Early','Late'},'ylim',[-6 26],'xticklabelrotation',45)
set(gca,'xTick',xticks,'xticklabel',barVarNames,'ylim',[-30 60],'xticklabelrotation',45);

%% BAR GRAPHS avg both directions
figure

% Plot bars
barVarNames = summaryBar.Properties.VariableNames;
for vi = 1:6
    dataPoints = summaryBar.(barVarNames{vi});
    dpPlotBar(vi, dataPoints );
end

xticks = [1:6];
set(gca,'xTick',xticks,'xticklabel',barVarNames,'ylim',[-30 60],'xticklabelrotation',45);

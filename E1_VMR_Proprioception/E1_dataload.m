tic
clear all; close all;
format compact

baseDir='C:\BoxSync\Dropbox\VICE\Darius\Prop_EC_pilot\Prop_EC_pilot\';
dataDir=fullfile(baseDir,'Data\E2');
analyzeDir=fullfile(baseDir,'Analyze\');

% List of participants
% Day 1 data first (name ends with A), then Day 2 (name ends with B)
subj = {'PropVMR_602A__S53A', 'PropVMR_602A__S54A', ... % begin day 1
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


T = []; % create data frame
M=[];

nt = 470; % number of trials

mm2pixel = 3.6137;  pixel2mm = 1/mm2pixel; pixel2cm = pixel2mm./10;
dt_tablet = 0.005;  % check to make sure this is correct
maxReachTime = 350; %700 ms, assuming 500 hz sampling rate

% Define some anonymous helper functions
NaNmask_ = @(x) (double(x)./double(x).*double(x));
Outlier_ = @(x,Niqr) (abs(x-nanmean(nanmedian(x))) > Niqr*nanmean(iqr(x)));


% for s = 1
for s = 1:length(subj) % loop participants
    subj{s} % pick participant
    load(fullfile(dataDir,[subj{s}, '.mat'])); % load participant's file
    
    
    %%%%% ANALYZE TRAJECTORIES BEGIN   
    
    S = [];
    
    S.hand_x = nan(nt,maxReachTime);     % game loop ran at 1000 hz
    S.hand_y = nan(nt,maxReachTime);
    S.hand_dist = nan(nt,maxReachTime);
    S.hand_v = nan(nt,maxReachTime);
    S.trialtime = nan(nt,maxReachTime);
    S.SN(1:nt,1) = s;
    S.tgtpos(1:nt,1) = [tgtloc(1,1)-xCenter];
    S.tgtpos(1:nt,2) = [yCenter-tgtloc(1,2)];
    
    V = [];
    Z = [];

    Z.move_trial = trial_move;
    Z.gamephase = gamephase_move;
       
    
    for i=1:nt
        
        hand_dist = sqrt(hand_x.^2 + hand_y.^2);
        
        % trimming hand data
        idx1 = find(Z.move_trial==i & (Z.gamephase==2 | Z.gamephase==3 | Z.gamephase==4)); %look for reaching during outbound traj, endpoint, and btwn blocks phases
        idx2 = find(Z.move_trial==i+1 & Z.gamephase==0,100); %participant may be reaching during start of next trial
        idx = [idx1;idx2];
        if length(idx)<=maxReachTime
            S.hand_x(i,1:length(idx)) = hand_x(idx);
            S.hand_y(i,1:length(idx)) = hand_y(idx);
            S.hand_dist(i,1:length(idx)) = hand_dist(idx);
            S.trialtime(i,1:length(idx)) = trial_time(idx);
        else
            S.hand_x(i,:) = hand_x(idx(1:maxReachTime));
            S.hand_y(i,:) = hand_y(idx(1:maxReachTime));
            S.hand_dist(i,:) = hand_dist(idx(1:maxReachTime));
            S.trialtime(i,:) = trial_time(idx(1:maxReachTime));
        end
        
    end
    
    [S.hx,S.hy,S.hdist,S.trialt] = deal(NaN(size(S.hand_x)));
    
    for k=1:nt
        
        % resample hand position based on samples where position
        % changes(this is ok because we're looking at data during movement
        % only)
        moveidx1 = abs(diff(S.hand_y(k,:))) > 0;  % creating logical to keep track of when there is movement
        moveidx2 = abs(diff(S.hand_x(k,:))) > 0;
        
        ii = [1, find(moveidx1+moveidx2)+1];    % indices of when there is new sample
        Nii(k) = length(ii);
        S.hx(k,1:Nii(k)) = S.hand_x(k,ii)*pixel2mm;
        S.hy(k,1:Nii(k)) = S.hand_y(k,ii)*pixel2mm + 80; % '+80' to compensate for the start position being at -80
        S.hdist(k,1:Nii(k)) = S.hand_dist(k,ii)*pixel2mm;
        S.trialt(k,1:Nii(k)) = S.trialtime(k,ii);
        
    end
    
    hvx = []; hvy = [];
    hvx = diff(S.hx')'./dt_tablet;    % compute velocity based on a simple difference (this can be smoothed later if necessary)
    hvy = diff(S.hy')'./dt_tablet;
    S.absvel = sqrt(hvx.^2 + hvy.^2);
    S.radvel = diff(S.hdist')'/dt_tablet;
    S.absacc = diff(S.absvel')'/dt_tablet';     % compute acceleration
    
    % Remove small number of huge velocity spikes -- figure this out!
    zz = zeros(1,nt);
    q = Outlier_([zz; diff(S.hdist')],10)';
    S.hdist_SpikeNum = sum(sum(q));
    S.hdist = S.hdist.*NaNmask_(~q);
    
    q = Outlier_([zz; diff(S.absvel')],10)';
    S.absvel_SpikeNum = sum(sum(q));
    S.absvel = S.absvel.*NaNmask_(~q);
    
    q = Outlier_([zz; diff(S.radvel')],10)';
    S.radvel_SpikeNum = sum(sum(q));
    S.radvel = S.radvel.*NaNmask_(~q);
    
    q = Outlier_([zz; diff(S.absacc')],10)';
    S.absacc_SpikeNum = sum(sum(q));
    S.absacc = S.absacc.*NaNmask_(~q);
    
    
    S.xi = S.hx(:,1);   % save some basic info about the movement: initial position, MT, max y-vel
    S.yi = S.hy(:,1);
    S.MT = (Nii/dt_tablet)';    
    [S.absvelmax, S.absvelmax_idx] = max(S.absvel');
    S.absvelmax = S.absvelmax';
    S.absvelmax_idx = S.absvelmax_idx';
    
    [S.radvelmax, S.radvelmax_idx] = max(S.radvel');
    S.radvelmax = S.radvelmax';
    S.radvelmax_idx = S.radvelmax_idx';
       
    for k = 1:nt
        
        S.xf(k,1) = S.hx(k,Nii(k));
        S.yf(k,1) = S.hy(k,Nii(k));
        
        S.maxv_hand_ang(k,1) = atan2d(S.hy(k,S.absvelmax_idx(k)), S.hx(k,S.absvelmax_idx(k)));
        
        S.maxradv_hand_ang(k,1) = atan2d(S.hy(k,S.radvelmax_idx(k)), S.hx(k,S.radvelmax_idx(k)));
        
        % Early hand angle(50 ms); based on sampling rate, 5 ms in between
        % samples
        S.hand_ang_50(k,1) = atan2d(S.hy(k,11), S.hx(k,11));
        
        actualtime_50(k,1) = S.trialt(k,11)-S.trialt(k,1);
        
        if sum(S.hdist(k,1:end-1)>=80 & S.radvel(k,:)<=0) > 0
            S.turnIdx(k,1) = find(S.hdist(k,1:end-1)>=80 & S.radvel(k,:)<=0,1);
            S.maxRadDist(k,1) = S.hdist(k,S.turnIdx(k,1));
            if S.turnIdx(k,1)-1 > 0  % DP: not sure why but one subject lacked data for this trial. used this to prevent it crashing when S.maxRadDist2(k,1) = 0
                S.maxRadDist2(k,1) = S.hdist(k,S.turnIdx(k,1)-1);
            else
                S.maxRadDist2(k,1)= nan;
            end
            S.handAngMaxDist(k,1) = atan2d(S.hy(k,S.turnIdx(k)), S.hx(k,S.turnIdx(k)));
        else
            S.turnIdx(k,1) = 0;
            S.maxRadDist(k,1) = nan;
            S.maxRadDist2(k,1) = nan;
            S.handAngMaxDist(k,1) = nan;
        end
        
        S.testMaxRadDist(k,1) = max(S.hdist(k,:));
        
    end
    
    % rotate target to zero to get hand angle
    for i = 1:nt
        
        theta(i) = atan2d(sind(hand_angle(i) - tgt_ang(i)), cosd(hand_angle(i) - tgt_ang(i)));
        
        % calculate hand angle at max velocity
        theta_maxv(i) = atan2d(sind(S.maxv_hand_ang(i) - tgt_ang(i)), cosd(S.maxv_hand_ang(i) - tgt_ang(i)));
        
        % calculate hand angle at peak radial velocity
        theta_maxradv(i) = atan2d(sind(S.maxradv_hand_ang(i) - tgt_ang(i)), cosd(S.maxradv_hand_ang(i) - tgt_ang(i)));
        
        % calculate hand angle at 100 ms after movement
        theta_50(i) = atan2d(sind(S.hand_ang_50(i) - tgt_ang(i)), cosd(S.hand_ang_50(i) - tgt_ang(i)));
        
        % calculate hand angle at maximum radial extent of reach
        thetaMaxExtent(i) = atan2d(sind(S.handAngMaxDist(i) - tgt_ang(i)), cosd(S.handAngMaxDist(i) - tgt_ang(i)));
    end
    
    % Compile relevant variables
    V.hx = S.hx;
    V.hy = S.hy;
    V.absvel = S.absvel;
    V.absacc = S.absacc;
    V.hdist = S.hdist;
    V.radvel = S.radvel;
    V.maxRadDist = S.maxRadDist;
    
    
    %%%%% ANALYZE TRAJECTORIES END
    
    
    %%%%% ANALYZE TRIALS BEGIN
    
    % this gets us hand angle after target rotated to zero
    for i = 1:nt
        theta(i) = atan2d(sind(hand_angle(i)-tgt_ang(i)), cosd(hand_angle(i)-tgt_ang(i)));
    end
    
    % save file number
    D.subj_filenum(1:nt,1) = str2num(subj{s}(16:17));
    
    
    if subj{s}(length(subj{s})) == 'A' % check if last letter in subject name is A (A is for Day 1)
        D.SN(1:nt,1) = s;% subject number, by index
        D.day_num(1:nt,1) = 1; % day 1
    elseif subj{s}(length(subj{s})) == 'B' % check if last letter in subject name is B (B is for Day 2)
        D.SN(1:nt,1) = T.SN(T.subj_filenum == nanmean(D.subj_filenum) & T.day_num == 1); % save same subject number as the file for day 1
        D.day_num(1:nt,1) = 2; % day 2
    end
    
    
    D.rot_cond(1:nt,1) = rotation(196); % rotation condition 
    
    D.prop_homex(1:nt,1) = (homex(1)-xCenter)*(1/mm2pixel); % x-val of home, in mm2 & with origin changed to (0,0)
    D.prop_homey(1:nt,1) = (homey(1)-yCenter)*(1/mm2pixel)*(-1); % y-val of home, in mm2 & with origin changed to (0,0)
    
        
    D.TN(1:nt,1) = (1:nt); % trial number
    D.hand_theta = theta'; 
    D.hand_theta_maxv = theta_maxv';
    D.hand_theta_maxradv = theta_maxradv';
    D.handMaxRadExt = thetaMaxExtent';
    D.hand_theta_50 = theta_50';
    D.raw_ep_hand_ang = hand_angle;
    D.ti = tgt_ang; % target angle
    D.fbi = online_fb; % feedback on or off
    D.ri = rotation; % rotation on or off
    D.MT = MTs'; % movement time
    D.RT = RTs'; % reaction time
    
    D.robot_TT = robotTimes; % movement time of researcher during prop trials (from home to target)
    
    D.PropTestAng = prop_test_ang; % prop target angle
    
    % NEW VERSION OF PROPLOCX
    D.PropLocX = (prop_testloc(:,1)-homex).*(1/mm2pixel); % x-val of prop test location, in mm2 & origin at (0,0)
    D.PropLocY = (prop_testloc(:,2)-homey)*(-1)*(1/mm2pixel); % Flipping y values because game has (0,0) at top left of screen

%     % OLD VERSION OF PROPLOCK - ALISSA HAD IT LIKE THIS. 
%     D.PropLocX_old = (prop_testloc(:,1)-xCenter).*(1/mm2pixel); % x-val of prop test location, in mm2 & origin at (0,0)
%     D.PropLocY_old = (prop_testloc(:,2)-yCenter)*(-1)*(1/mm2pixel); % Flipping y values because game has (0,0) at top left of screen
    
    
    D.PB = prop_block; % prop block
    D.RB = nan(nt,1); % reaching block (filled in later)
   
    
    D.FC_TT = FCTimes; % trial time for prop
    
    D.FC_X = (clickPoints(:,1)-homex)*(1/mm2pixel); % x-val of prop choice, in mm2 & origin at (0,0)
    D.FC_Y = (clickPoints(:,2)-homey)*(-1)*(1/mm2pixel); % Flipping y values because game has (0,0) at top left of screen 
    D.HL_X = (hL_X-homex)*(1/mm2pixel); % hand location during click, in mm2 & origin at (0,0)
    D.HL_Y = (hL_Y-homey)*(-1)*(1/mm2pixel); % Flipping y values because game has (0,0) at top left of screen
        
    D.FC_bias_X = (clickPoints(:,1)-hL_X)*(1/mm2pixel); % distance between actual hand location & guessed location, in mm2 & origin at (0,0) 
    D.FC_bias_Y = (clickPoints(:,2)-hL_Y)*(1/mm2pixel)*(-1); % Flipping y values because game has (0,0) at top left of screen
    
    D.FC_X = (clickPoints(:,1)-homex)*(1/mm2pixel); % x-val of prop choice, in mm2 & origin at (0,0)
    D.FC_Y = (clickPoints(:,2)-homey)*(-1)*(1/mm2pixel); % Flipping y values because game has (0,0) at top left of screen 
    D.HL_X = (hL_X-homex)*(1/mm2pixel); % hand location during click, in mm2 & origin at (0,0)
    D.HL_Y = (hL_Y-homey)*(-1)*(1/mm2pixel); % Flipping y values because game has (0,0) at top left of screen
        
    D.FC_bias_X = (clickPoints(:,1)-hL_X)*(1/mm2pixel); % distance between actual hand location & guessed location, in mm2 & origin at (0,0) 
    D.FC_bias_Y = (clickPoints(:,2)-hL_Y)*(1/mm2pixel)*(-1); % Flipping y values because game has (0,0) at top left of screen
    
    % prop_theta this gets us hand angle after target rotated to zero
    for i = 1:nt
        HL_angle(i,1) = atan2d( D.HL_Y(i,1), D.HL_X(i,1) ); % Aug 28 2018: "-80" is so that the angle is worked out relative to the actual start position
        FC_angle(i,1) = atan2d( D.FC_Y(i,1), D.FC_X(i,1) ); % Aug 28 2018: "-80" is so that the angle is worked out relative to the actual start position
        
        prop_theta(i,1) = atan2d(sind(FC_angle(i)-HL_angle(i)), cosd(FC_angle(i)-HL_angle(i)));
    end
    D.prop_theta = prop_theta;
    
    %%%%% ANALYZE TRIALS END
    
  
%     % Remove the clamp demonstration practise trials (only need to remove these for error
%     clamp experiments)
%     D = structfun(@(x) (x([(1:trialbeforepractice),(startofclamp:end)],:)), D, 'uniformoutput', 0);    
%     V = structfun(@(x) (x([(1:trialbeforepractice),(startofclamp:end)],:)), V, 'uniformoutput', 0);
    
    % Reassign TN to account for removed practise trials
    D.TN(1:length(D.TN),1) = 1:length(D.TN);
    D.CN(1:length(D.TN),1) = ceil((1:length(D.TN))/5)';
    
    % Create RB for cycle number
    reaching_blockCN = [1 2; 3 4; 15 32; 33 33; 40 47; 54 61; 62 62; 69 76; 77 77; 84 84 ; 85 94];
    for bi = 1:length(reaching_blockCN)
        CN_idx = (D.CN>=reaching_blockCN(bi,1) & D.CN<=(reaching_blockCN(bi,2)));
        D.RB(CN_idx) = bi  ;
    end
    
    % Cleaning up some variables (should be nan rather than 0)
    D.ti(D.ti==0) = nan; D.PropTestAng(D.PropTestAng==0) = nan; D.PB(D.PB==0) = nan; 
    
    
    
    % Compile final tables
    
    temp = struct2table(D);
    tempmove = struct2table(V);
    
    T = [T;temp];
    M = [M;tempmove];
    
    elapsed_times(s) = elapsedTime;
end

mean_elapsed_experiment_time = mean(elapsed_times)/60;
minutes = floor(mean_elapsed_experiment_time)
seconds = 60*(mean_elapsed_experiment_time  - floor(mean_elapsed_experiment_time) )

cd(analyzeDir)
save('E1_trials.mat','T');
save('E1_movefile.mat','M');
toc

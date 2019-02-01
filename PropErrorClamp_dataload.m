% This function loads in all subject data.

tic
clear all; close all

baseDir='C:\BoxSync\Dropbox\VICE\Darius\PropErrorClamp\';
dataDir=fullfile(baseDir,'Data\');
analyzeDir=fullfile(baseDir,'Analyze\');

subj = {'PropEC_622C__J9','PropEC_622A__ j5_a',... % First two are RAs that were usable data
        'PropEC_905C_1_J1', 'PropEC_905A_1_J2_a', 'PropEC_905C__J3_a', 'PropEC_905A__J4_a',...
        'PropEC_905C__J5_a', 'PropEC_905A__J6_a', 'PropEC_905C__J7_a', 'PropEC_905A__J8_a',...
        'PropEC_905C__J9_a','PropEC_905A__J10_a', 'PropEC_905C__J11_a','PropEC_905A__J12_b',...
        'PropEC_905C__J13_a', 'PropEC_905A__J14', 'PropEC_905C__J15_a', 'PropEC_905A__J16_a',...
        'PropEC_905C__J17_a', 'PropEC_905A__J18_a', 'PropEC_905C__J19_c','PropEC_905A__J20',...
        'PropEC_905C__J21_a','PropEC_905A__J22_a','PropEC_905C__J23','PropEC_905A__J24_a',...
        'PropEC_905C__J25_a','PropEC_905A__J18_b','PropEC_905A__J26_a', 'PropEC_905C__J27_a',...
        'PropEC_905C__J28_a','PropEC_905A__J29'};
    
%%%%% Practise trials get removed later %%%%%
trialbeforepractice = 180;
startofclamp = 184;


%%% look for table. if it's there, we will append stuff. look to see how
%%% many subjects are in existing table. 
if (exist('pre-processed_data\PropErrorClamp_trials.mat') & exist('pre-processed_data\PropErrorClamp_movefile.mat')) == 1
    load('pre-processed_data\PropErrorClamp_trials.mat'); load('pre-processed_data\PropErrorClamp_movefile.mat')
    T=T;
    M=M;
    num_tested_subj=numel(unique(T.SN));
    if num_tested_subj==length(subj)
        return
    end
else
    T=[];
    M=[];
    num_tested_subj=0;
end

T=[];
M=[];

nt = 777;   % number of trials
tpe = 9;    % trials per epoch  CHECK HOW MANY TRIALS PER BLOCK

mm2pixel = 3.6137;  pixel2mm = 1/mm2pixel; pixel2cm = pixel2mm./10;
dt_tablet = 0.005;  % check to make sure this is correct
maxReachTime = 350; %700 ms, assuming 500 hz sampling rate

% Define some anonymous helper functions
NaNmask_ = @(x) (double(x)./double(x).*double(x));
Outlier_ = @(x,Niqr) (abs(x-nanmean(nanmedian(x))) > Niqr*nanmean(iqr(x)));

% loop through all subjects
for s = num_tested_subj+1:length(subj)
        
    s
    subj{s}
    load(fullfile(dataDir,[subj{s}, '.mat']));
    
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
        S.hy(k,1:Nii(k)) = S.hand_y(k,ii)*pixel2mm;
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
    
    
    % may need to change these lines when testing zero clamp or change in clamp sizes   
    if max(clamped_feedback)==1 & max(rotation) > 0
        direction = 1;
        condition = max(rotation);
    elseif max(clamped_feedback)==1 & max(rotation) <= 0
        direction = 0;
        condition = min(rotation);
    end
    
    targetsize = targetsize*pixel2mm;
    grp=abs(condition);    
    
   
    
    
    D.SN(1:nt,1) = s;
    D.tester(1:nt,1) = subj{s}(1);
    D.group(1:nt,1) = grp;
    D.TN(1:nt,1) = nan(nt,1); % trial number
    D.CN(1:nt,1) = nan(nt,1); % cycle number
    D.BN(1:nt,1) = nan(nt,1); % block number
    D.rot_cond(1:nt,1) = condition;
    D.CCW(1:nt,1) = direction;
    D.tgtsize(1:nt,1) = targetsize;
    D.hand_theta = theta';
    D.hand_theta_maxv = theta_maxv';
    D.hand_theta_maxradv = theta_maxradv';
    D.handMaxRadExt = thetaMaxExtent';
    D.hand_theta_50 = theta_50';
    D.raw_ep_hand_ang = hand_angle;
    D.ti = tgt_ang;
    D.fbi = online_fb;
    D.ri = rotation;
    D.clampi = clamped_feedback;
    D.MT = MTs';
    D.RT = RTs';
    D.ST = SearchTimes';
    D.radvelmax = S.radvelmax;
    D.maxRadDist = S.maxRadDist;
    D.testMaxRadDist = S.testMaxRadDist;  
    
    % Variables specific to Proprioceptive Error Clamp experiment
    D.robot_TT = robotTimes; % movement time of researcher during prop trials (from home to target)
        
    D.PropTestAng = prop_test_ang; % prop target angle
    D.PropLocX = (prop_testloc(:,1)-homex).*(1/mm2pixel); % x-val of prop test location, in mm2 & origin at (0,0)
    D.PropLocY = (prop_testloc(:,2)-homey)*(-1)*(1/mm2pixel); % Flipping y values because game has (0,0) at top left of screen

    D.PB = prop_block; % prop block
    D.RB = nan(nt,1); % reaching block (filled in later)
        
    D.FC_TT = FCTimes; % trial time for prop
    
    % Variables relating to estimated hand position and proprioceptive bias
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
    
    
    % Trajectory information  
    V.hx = S.hx;
    V.hy = S.hy;
    V.absvel = S.absvel;
    V.absacc = S.absacc;
    V.hdist = S.hdist;
    V.radvel = S.radvel;
    V.maxRadDist = S.maxRadDist;
    
    
    % Remove the practise trials
    D = structfun(@(x) (x([(1:trialbeforepractice),(startofclamp:end)],:)), D, 'uniformoutput', 0);    
    V = structfun(@(x) (x([(1:trialbeforepractice),(startofclamp:end)],:)), V, 'uniformoutput', 0);
    
    % Reassign TN and move_cycle to account for removed practise trials   
    D.TN(1:length(D.TN)) = 1:length(D.TN);
    
    % Create CN for cycle number
    D.CN(1:length(D.TN)) = ceil((1:length(D.TN))/tpe)';
    
    % Create reaching block number (RB) for cycle number
    reaching_blockCN = [2 12; 20 40; 44 54; 58 68; 72 82];
    for bi = 1:length(reaching_blockCN)
        CN_idx = (D.CN>=reaching_blockCN(bi,1) & D.CN<=(reaching_blockCN(bi,2)));
        D.RB(CN_idx) = bi  ;
    end
    
    
    % Cleaning up some variables (should be nan rather than 0)
    D.hand_theta_maxv(D.ti==0) = nan; D.hand_theta_maxradv(D.ti==0) = nan;
    D.handMaxRadExt(D.ti==0) = nan; D.hand_theta_50(D.ti==0) = nan;
    D.ti(D.ti==0) = nan; D.PropTestAng(D.PropTestAng==0) = nan; D.PB(D.PB==0) = nan; 
    D.hand_theta_maxv(D.ti==0) = nan;
    
    
    temp = struct2table(D);
    tempmove = struct2table(V);
    
    T = [T;temp];
    M = [M;tempmove];
    
    elapsed_times(s) = elapsedTime;
    
end

mean_elapsed_experiment_time = mean(elapsed_times)/60;
minutes = floor(mean_elapsed_experiment_time)
seconds = 60*(mean_elapsed_experiment_time  - floor(mean_elapsed_experiment_time) )


save('pre-processed_data\PropErrorClamp_trials.mat','T');
save('pre-processed_data\PropErrorClamp_movefile.mat','M');
toc


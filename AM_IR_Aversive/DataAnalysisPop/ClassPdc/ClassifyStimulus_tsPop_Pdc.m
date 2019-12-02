function ClassifyStimulus_tsPop
% ClassifyStimulus_SU
%
%  Based on RCorr_SU. Changing decision rule to euclidean distance, which
%  is applied in rc_CalcEucDist [was rc_CalcR]. Middle-man function is still
%  run_RCorr.
%  
%
%  KP, 2019-10
%


close all

TempSize    = 1;
app_str     = '_1trTemp';
StimDur     = 500;

convwin     = 10;
convwin     = ones(1,convwin).*(1/convwin);


%% Load data

fn = set_paths_directories('','',1);

% Load spikes data (created in cumulativeSpikeCount)
load(fullfile(fn.figs,'CmSpCount',['Cell_Time_Trial_Stim_' num2str(20) 'trs']))

% Load Unit data files
fn = set_paths_directories('','',1);
q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q

% % % Find just Apr02 cells
% % sess_idx = find(strcmp(UnitInfo.Subject,'AAB_265054') & strcmp(UnitInfo.Session,'Apr02-AM'));
% % Cell_Time_Trial_Stim = Cell_Time_Trial_Stim(intersect(Un_Indices,sess_idx),:,:,:);


%% Preallocate

Cell_Time_Trial_Stim = Cell_Time_Trial_Stim(:,1:StimDur,:,:);
nUnits  = size(Cell_Time_Trial_Stim,1);
nTrials = size(Cell_Time_Trial_Stim,3);
nStim   = size(Cell_Time_Trial_Stim,4);

PCMat = nan(nStim,nStim);
percentCorrectMat = nan(nStim,nStim,nTrials);

% figure; hold on

for iIt = 1:nTrials
    
    % Set templates
    T = cell(2,nStim);
    for ist = 1:nStim
        
        T(1,ist) = {sum(Cell_Time_Trial_Stim(:,:,iIt,ist),1)};
%         T(1,ist) = {reshape(Cell_Time_Trial_Stim(:,:,iIt,ist),1,nUnits*StimDur,1,1)};
        foo = conv(T{1,ist},convwin);
        T(1,ist) = {foo(floor(numel(convwin)/2)+(0:StimDur-1))};
%         subplot(3,2,ist); hold on
%         plot(T{1,ist})
%         T(1,ist) = {Cell_Time_Trial_Stim(:,:,iIt,ist)'};
    end
    T(2,:) = deal({iIt});
    
    % Go through each stimulus
    for isTrue = 1:nStim
        
        % Preallocations
        finalMat = [];
        count = 0;
        
        % Go through each trial
        test_trials = 1:nTrials;
        test_trials(iIt) = [];
        for it = test_trials
            
            S = sum(Cell_Time_Trial_Stim(:,:,it,isTrue),1);
%             S = reshape(Cell_Time_Trial_Stim(:,:,it,isTrue),1,nUnits*StimDur,1,1);
            foo = conv(S,convwin);
            S = foo(floor(numel(convwin)/2)+(0:StimDur-1));
            
            [blockAssignment,maxR] = rc_calcEucDist(T,S); %min Euc dist
%             [blockAssignment,maxR] = rc_calcR(T,S);      %max rcorr
            count = count+1;
            finalMat(count,:) = [it blockAssignment maxR];
            
        end %it
        
        for isAss =  1:nStim
            percentCorrectMat(isTrue,isAss,iIt) = sum(finalMat(:,2)==isAss)/size(finalMat,1);
        end
        
    end %isTrue
    
end %iIt

PCMat = mean(percentCorrectMat,3,'omitnan');



%% Save data

savedir = fullfile(fn.figs,'StimClass');
% if exclOnset>0
%     savedir = fullfile(savedir,'exclOnset');
% end
if ~exist(savedir,'dir')
    mkdir(savedir)
end

savename = ['PCMat_tsPopAvg_conv10_Ed' app_str '_' num2str(StimDur) 'ms'];

save(fullfile(savedir,savename),'PCMat','-v7.3')


%% Plot results

% Figure settings
set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',18)

% scrsz = get(0,'ScreenSize');     %[left bottom width height]
% fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];
% shortscreen = [1 scrsz(4)/3 scrsz(3) scrsz(4)/2];

hf2 = figure;
imagesc(PCMat)
axis square
caxis([0 1])
cmocean('ice') %ice
% cmocean('curl','pivot',0)
colorbar
ylabel('True stim')
xlabel('Assigned')

% print_eps_kp(hf2,fullfile(savedir,['PCMat_nSpk_Pop' app_str '_' num2str(StimDur) 'ms']))
print(fullfile(savedir,savename),'-dpdf')


end


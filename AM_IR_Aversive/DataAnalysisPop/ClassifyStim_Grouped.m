function ClassifyStim_Grouped
% ClassifyStim_Grouped
% 
%  Linear classifier for population decoding, concatenated summed activity 
%  within groups of cells.
% 
%  KP, 2019-11
%

% 1 trial
%
%  RESULTS ARE ROUGHLY THE SAME RANGE, NO MATTER THE GROUPING TYPE
%  OR THE DECISION RULE. SMOOTHING THE DATA IMPROVES PC A BIT. 

%  It seems like the single trial linear euclidean classifier is not
%  sufficient to discern differences between grouping methods, because
%  running too close to the floor. Could try with larger templates, to see
%  if it's worth using a fancier classifier.

%  A single group (3) does better than all concatenated. 


% 10 trials
%
%  Also similar performance. 
% 


%~~~~~~~~~~~~~~~~~~~~~~~~~~
binsize   =  10;
%~~~~~~~~~~~~~~~~~~~~~~~~~~
nGrps     =  5;
GrpType   =  'dynRange'; 'tuning'; 'peakFR'; 'phase'; 
Randomize =  0;
%~~~~~~~~~~~~~~~~~~~~~~~~~~
SetNIterations = 500;
%~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_ex   =  0;
%~~~~~~~~~~~~~~~~~~~~~~~~~~
TempSize  =  19;
app_str   =  '_19trTemp';
%~~~~~~~~~~~~~~~~~~~~~~~~~~
StimDur   =  500;
%~~~~~~~~~~~~~~~~~~~~~~~~~~

convwin   =  ones(1,binsize).*(1/binsize);
rng('shuffle')


%% Load Unit data files

fn = set_paths_directories('','',1);

q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q
%-------
spkshift = mean([UnitData([UnitData.IntTime_spk]>0).IntTime_spk]);
%-------

% Load specified groupings of cells
grpfn = sprintf('GroupLabels_%s_%i.mat',GrpType,nGrps);
load(fullfile(fn.figs,'PopMPH','Grouped',grpfn))


% Get spiking data 
% First, easy way: use the one set of 20 random trials already saved
% Load spikes data (created in cumulativeSpikeCount)
load(fullfile(fn.figs,'CmSpCount',['Cell_Time_Trial_Stim_' num2str(20) 'trs']))

if size(GroupLabel,1) ~= max(Un_Indices)
    keyboard
end


%% Figure settings

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',18)

scrsz = get(0,'ScreenSize');     %[left bottom width height]
widescreen  = [1 scrsz(4)/2 scrsz(3) scrsz(4)/2];
% shortscreen = [1 scrsz(4)/3 scrsz(3) scrsz(4)/2];

% Set colors
colors = [ 150 150 150;...
    84  24  69;...
    120  10  41;...
    181   0  52;...
    255  87  51;...
    255 153   0]./255;

savedir = fullfile(fn.figs,'StimClass','GroupEach');
if ~exist(savedir,'dir')
    mkdir(savedir)
end


%% Group data for input

nTrials  = size(Cell_Time_Trial_Stim,3);
nStim    = size(Cell_Time_Trial_Stim,4);

Cell_Time_Trial_Stim  = Cell_Time_Trial_Stim(:,1:StimDur,:,:);

Group_Time_Trial_Stim  = nan( max(GroupLabel), size(Cell_Time_Trial_Stim,2), ...
                         size(Cell_Time_Trial_Stim,3), size(Cell_Time_Trial_Stim,4) );


%% 
% RandGrpPCMat = nan(nStim,nStim,1000);
% 
% for iii = 1:size(RandGrpPCMat,3)
                     
allUnits   = 1:numel(GroupLabel);
allUnsUsed = [];

for iGrp = 1:max(GroupLabel)
    GrpUns = find(GroupLabel==iGrp);
    
    % Randomize, if desired
    if Randomize
        GrpUns = allUnits(randperm(numel(allUnits),numel(GrpUns)));
        allUnits(ismember(allUnits,GrpUns)) = [];
        allUnsUsed = [allUnsUsed GrpUns];
        GrpType   = 'random';
    end
    
    iCTTS = ismember(Un_Indices,GrpUns);
    Group_Time_Trial_Stim(iGrp,:,:,:) = sum(Cell_Time_Trial_Stim(iCTTS,:,:,:),1);
end


%% Classify

nGroups  = size(Group_Time_Trial_Stim,1);

if TempSize==1 || TempSize==(nTrials-1)
    nIterations = nTrials;
else
    nIterations = SetNIterations;  %new random draws of trials for templates
end

if plot_ex
    
    figure; hold on
% %     set(gcf,'Position',widescreen)
% %     
% %     xplot = 1:500;
% %     
% %     for ist = 1:nStim
% % %         subplot(3,2,ist);
% % %         hold on
% %         for ig = 1:size(Group_Time_Trial_Stim,1)
% %             data = mean(Group_Time_Trial_Stim(ig,:,:,ist),3);
% %             foo  = conv(data,convwin);
% %             data = foo(floor(numel(convwin)/2)+(0:StimDur-1));
% %             plot(xplot+max(xplot)*(ig-1),data,'LineWidth',2,'Color',colors(ist,:))
% %         end
% %     end
% %     title(['Cells grouped by ' GrpType])
% %     set(gca,'Color','none','tickdir','out',...
% %         'xtick',0:250:2500,'xticklabel',{'' 'Group 1' '' 'Group 2' '' 'Group 3' '' 'Group 4' '' 'Group 5' '' })
% %     
% %     % Save figure
% %     savename = sprintf('ExResponses_grp_%s5_%ims', GrpType,binsize);
% %     print_eps_kp(gcf,fullfile(savedir,savename))
    
end


percentCorrectMat = nan(nStim,nStim,nIterations);

for iIt = 1:nIterations
    
    % Choose template and test trials
    if TempSize==1
        temp_trs = iIt;
    elseif TempSize==(nTrials-1)
        temp_trs = 1:nTrials; temp_trs(iIt)=[];
    else
        temp_trs = randperm(nTrials,TempSize);
    end
    test_trs = 1:nTrials; test_trs(temp_trs)=[];
    
    
    % Set templates
    T = cell(2,nStim);
    
    for ist = 1:nStim
        
        % Concatenate data from each group
        T(1,ist) = {reshape(mean(Group_Time_Trial_Stim(:,:,temp_trs,ist),3),1,nGroups*StimDur,1,1)};
        
        % Bin to something larger than 1 ms
        %         T(1,ist) = {mean(reshape(T{1,ist},binsize,size(T{1,ist},2)./binsize),1)};
        
        % Instead of binning, convolve with boxcar
        foo = conv(T{1,ist},convwin);
        T(1,ist) = {foo(floor(numel(convwin)/2)+(0:StimDur-1))};
        
        
        if plot_ex
            subplot(3,2,ist); hold on
            plot(T{1,ist})
            xlim([1 size(T{1,ist},2)])
            ylim([0 2])
        end
        
    end
    T(2,:) = deal({iIt});
    
    
    % Go through each stimulus
    for isTrue = 1:nStim
        
        % Preallocations
        finalMat = [];
        count = 0;
        
        % Step through each trial that is not part of the template
        
        for it = test_trs
            
            % Concatenate data from each group
            S = reshape(Group_Time_Trial_Stim(:,:,it,isTrue),1,nGroups*StimDur,1,1);
            
            % Bin to something larger than 1 ms
            %             S = mean(reshape(S,binsize,size(S,2)./binsize),1);
            
            % Instead of binning, convolve with boxcar
            foo = conv(S,convwin);
            S = foo(floor(numel(convwin)/2)+(0:StimDur-1));
            
            
            [blockAssignment,maxR] = rc_calcEucDist(T,S); %min Euc dist
%             [blockAssignment,maxR] = rc_calcR(T,S);      %max rcorr
            
            if plot_ex
                subplot(3,2,blockAssignment); hold on
                plot(S,'r')
            end
            
            count = count+1;
            finalMat(count,:) = [it blockAssignment maxR];
            
        end %it
        
        for isAss =  1:nStim
            percentCorrectMat(isTrue,isAss,iIt) = sum(finalMat(:,2)==isAss)/size(finalMat,1);
        end
        
    end %isTrue
    
end %iIt

PCMat = mean(percentCorrectMat,3,'omitnan');

% RandGrpPCMat(:,:,iii) = PCMat;
% end %iii
% 
% PCMat  = mean(RandGrpPCMat,3);


% Extract statistics from confusion matrix
meanPC = mean(diag(PCMat));
dprime = norminv(mean(diag(PCMat)),0,1) - norminv(mean(PCMat(~diag(diag(PCMat)))),0,1); 


%% Plot results
    
hf2 = figure;
imagesc(PCMat)
axis square
caxis([0 1])
cmocean('thermal') %ice
% cmocean('curl','pivot',0)
colorbar
ylabel('True stim')
xlabel('Assigned')
title(sprintf('%s %i  |  mean PC=%0.3f', GrpType,iGrp,meanPC)) 


%% Save data

% Save figure & results 
savename = sprintf('PCMat_grp_%s_%ims_tempsz%i',GrpType,binsize,TempSize);
print(fullfile(savedir,savename),'-dpdf')
save(fullfile(savedir,savename),'PCMat','-v7.3')


end


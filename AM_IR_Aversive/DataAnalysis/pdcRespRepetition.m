function pdcRespRepetition
%  
% KP, 2019-04
%

close all
minNtr = 10;

%% Load Unit data files

fn = set_paths_directories('','',1);

q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q
%-------
spkshift = mean([UnitData([UnitData.IntTime_spk]>0).IntTime_spk]);
%-------

AMrates = [2 4 8 16 32];

%% Figure settings

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',14)

scrsz = get(0,'ScreenSize');   %[left bottom width height]
tallfig  = [1 scrsz(4) scrsz(3)/2 scrsz(4)];

USE_MEASURE = 'FR';
switch USE_MEASURE
    case 'FR'
        units = 'Hz';
        switch units
            case 'Hz'
                axmin = 0.01;
                axmax = 10*ceil(max([UnitData.BaseFR])/10);
            case 'z'
                axmin = -1;
                axmax = 1.5;
        end
    case {'TrV' 'FF'}
        units = ' ';
        axmin = 0.01;
        axmax = 4;
end


% Set directory for saving figs 
savedir = fullfile(fn.figs,'PdcRepetition');
if ~exist(savedir,'dir')
    mkdir(savedir)
end


%% Select a datapoint (or a few?) to highlight

ex_subj = 'AAB_265058';
ex_sess = 'Jan17-AM';
ex_ch   = 63;
ex_clu  = 1230;

ex_Un   = find(strcmp({UnitData.Subject},ex_subj) & strcmp({UnitData.Session},ex_sess) & [UnitData.Channel]==ex_ch & [UnitData.Clu]==ex_clu);
ex_Un   = numel(UnitData);

%% 

for iUn = 1:numel(UnitData)
    
    close all
    
    % -> -> -> -> -> -> -> -> -> ->
    %%% skips merged units for now
    if numel(UnitInfo(iUn,:).Session{:})==4  %strncmp(UnitInfo.RespType{iUn},'merged',6)
        continue
    end
    
    subject     = UnitData(iUn).Subject;
    session     = UnitData(iUn).Session;
    channel     = UnitData(iUn).Channel(1);
    clu         = UnitData(iUn).Clu(1);
    
    % Get sound parameters
    dBSPL       = UnitData(iUn).spl;
    LP          = UnitData(iUn).lpn;
        
    
    % Load data files
    
    if (iUn>1 && ~( strcmp(subject,UnitData(iUn-1).Subject) && strcmp(session,UnitData(iUn-1).Session) )) || iUn==1 ||  ~exist('Phase0','var') 
        fprintf('Loading %s sess %s...\n',subject,session)
        clear Info TrialData SoundStream SpoutStream RateStream Phase0
        filename = sprintf( '%s_sess-%s_Info'     ,subject,session); load(fullfile(fn.processed,subject,filename));
        filename = sprintf( '%s_sess-%s_TrialData',subject,session); load(fullfile(fn.processed,subject,filename));
    end
    if (iUn>1 && ~( strcmp(subject,UnitData(iUn-1).Subject) && strcmp(session,UnitData(iUn-1).Session) && channel==UnitData(iUn-1).Channel ) )  || iUn==1 || ( ~exist('Spikes','var') || ~exist('Clusters','var') )
        clear Clusters Spikes
        filename = sprintf( '%s_sess-%s_Spikes'   ,subject,session); load(fullfile(fn.processed,subject,filename));
    end
    
    
    % Get spiketimes and shift based on calculated integration time 
    if exist('Spikes','var')                                 % >>> UMS <<<
        
        spiketimes = unique(round( Spikes.sorted(channel).spiketimes(Spikes.sorted(channel).assigns==clu') * 1000 + spkshift ));  %ms
        
    elseif exist('Clusters','var')                            % >>> KS <<<
        
        iClu = find([Clusters.maxChannel] == channel & [Clusters.clusterID] == clu);
        spiketimes = unique(round( Clusters(iClu).spikeTimes * 1000 + spkshift ))';
        
    end
    
    fprintf(' analyzing ch %i clu %i\n',channel,clu)
    
    
    
    %% Get data
    
    Slopes_FR = nan(1,numel(AMrates));
    Slopes_VS = nan(1,numel(AMrates));
    Slopes_Ph = nan(1,numel(AMrates));
    DeltaNspk = cell(1,numel(AMrates));
    
    for irate = AMrates
        
        % -> -> -> -> -> -> -> -> -> -> -> -> -> -> -> -> -> -> 
        % Right away, skip units with v low FR for this stimulus
        if nanmean(UnitData(iUn).FR_raw_tr(:,1+find(AMrates==irate))) < 2  %lower? iUn=155 4hz, eg
            continue
        end
        
        Rastermats = cell(irate,1);
        [Rastermats{:}] = deal(nan(0,ceil(1000/irate)));
        
        P0idx     = find(Phase0(2,:)==irate);
        Repeats   = P0idx(diff(P0idx)==1);
        PdOnsets  = Phase0(1,Repeats);
        
        ipd=1;
        while ipd <= numel(PdOnsets)
            
            Nreps = find(diff(Repeats(ipd:end))~=1,1,'first');
                        
            % for each period in this sequence
            for irep = 0:Nreps
                
                if (irep+1)>numel(Rastermats), continue, end
                
                ii = ipd+irep;
                
                % -> -> -> -> -> -> -> -> -> -> -> -> 
                if PdOnsets(ii)>TrialData.offset(end) 
                    break %cant verify sound params after last official trial (unless go back to raw data)
                end
                
                PdOffset  = floor(PdOnsets(ii) + 1000/irate)-1;
                TDidx     = find([TrialData.onset]<=PdOnsets(ii) & [TrialData.offset]>PdOnsets(ii));
                if isempty(TDidx)
                    TDidx = find([TrialData.onset]>=PdOnsets(ii),1,'first');
                end
                
                % Check that params steady and subject on spout and no artifact
                if  (   TrialData(TDidx,:).SPL == UnitData(iUn).spl     && ...
                        TrialData(TDidx,:).LP  == UnitData(iUn).lpn     && ...
                        ~ismember(TDidx, Info.artifact(channel).trials) && ...
                        mean(SpoutStream(PdOnsets(ii):PdOffset))>0.999    )
                    
                    % Add to corresponding raster
                    sp=[];
                    sp = spiketimes(spiketimes>=PdOnsets(ii) & spiketimes<=PdOffset) - PdOnsets(ii) +1;
                    rasterrow = zeros(1,ceil(1000/irate));
                    rasterrow(sp) = 1;
                    Rastermats{irep+1}(end+1,:) = rasterrow;
                    
                else
                    % -> -> -> -> -> -> -> -> -> ->
                    % continue to next first period
                    break
                end
            end
            
            % Skip ipd ahead to next first period
            ipd = ipd + Nreps;
            
        end %ipd
        
        
        % Quick check of n trials from beginning to end of periodic
        if (size(Rastermats{1},1) / size(Rastermats{end},1)) >1.5
            keyboard
        end
        
        
        %% ~~~~~~~~~~~~  MEASURE RESPONSES  ~~~~~~~~~~~~
        
        these_pds = find(cellfun(@(x) size(x,1),Rastermats)>minNtr)'; %lower?
        
        % -> -> -> -> -> -> ->
        if numel(these_pds)<2
            continue
        end
        
        % Set periods constituting beginning and end
        if irate==2
            x_beg = 1;
            x_end = 2;
        else
            x_beg  = 1:ceil(irate/4);         %first 250 ms
            x_end  = ceil(3*irate/4+1):irate; %last 250 ms
        end
        
        minTrs = min(cellfun(@(x) size(x,1),Rastermats));
        
        % Preallocate
        x   = [];  FR  = [];  
        VS  = [];  RS  = []; 
                   Ph  = [];  
                   
        for irep = 1:numel(these_pds)
            
            x = [x; irep*ones(size(Rastermats{irep},1),1)];
            
            %---- FR ----
            FR = [FR; sum(Rastermats{irep},2)*irate]; 
            
            %---- VS ----
            [tr,Spktime] = find(Rastermats{irep}); 
            [thisVS,thisRS,p] = vectorstrength(Spktime',1000/irate);
            if p<0.05 || thisRS>13.1
                VS = [VS; thisVS]; 
            else
                VS = [VS; nan]; 
            end
            RS = [RS; thisRS];
            
            %---- Mean phase ----
            mu = meanphase(Spktime',1000/irate);
            Ph  = [Ph; mu];
            
        end
        
        
        %~~~~~~~~ QUANTIFY EFFECT OF REPETITION IN THIS UNIT ~~~~~~~~
        
        % Linear fit of all datapoints
        c_FR = polyfit(x,FR,1);
        c_VS = polyfit(unique(x),VS,1);
        c_Ph = polyfit(unique(x),Ph,1);
        
        % Store slopes
        Slopes_FR(AMrates==irate) = c_FR(1);
        Slopes_VS(AMrates==irate) = c_VS(1);
        Slopes_Ph(AMrates==irate) = c_Ph(1);
        
        
        % Compare beginning and end
        Nspk_paired = nan(minTrs,2);
        % Beg
        nspks = [];
        nspks = cellfun(@(x) sum(x(1:minTrs,:),2), Rastermats(x_beg),'UniformOutput',false);
        Nspk_paired(:,1) = sum(horzcat(nspks{:}),2);
        % End
        nspks = [];
        nspks = cellfun(@(x) sum(x(1:minTrs,:),2), Rastermats(x_end),'UniformOutput',false);
        Nspk_paired(:,2) = sum(horzcat(nspks{:}),2);
        
        p_wsr = signrank(Nspk_paired(:,1),Nspk_paired(:,2));
%         p_kwt = kruskalwallis(Nspk_paired);
        p_rs  = ranksum(Nspk_paired(:,1),Nspk_paired(:,2));
        DeltaNspk{AMrates==irate} = [median(Nspk_paired,1) p_wsr];
        
        
        %% Individual unit plots 
        
        % (set threshold high to skip)
        
        %------------ FR ------------
        if iUn==ex_Un || abs(c_FR(1))>100
            
            figure; hold on
%             plotSpread(FR,'distributionIdx',x,'distributionColors','k','showMM',3)
%             hold on
%             plot(unique(x),polyval(c_FR,unique(x)),'r-','LineWidth',2)
            for irep = 1:irate
                fill(irep+[-0.2 0.2 0.2 -0.2],...
                    [quantile(FR(x==irep),0.25) quantile(FR(x==irep),0.25) quantile(FR(x==irep),0.75) quantile(FR(x==irep),0.75)],...
                    0.7.*[1 1 1],'EdgeColor','none')
                plot(irep,median(FR(x==irep)),'.b','MarkerSize',80)
            end
            
            ymaxval = 12*ceil( max(UnitData(iUn).FR_raw_tr(:,1+find(AMrates==irate))) /10 );
            xlim([0 irep+1])
            ylim([0 ymaxval])
            xlabel('# in periodic sequence')
            ylabel('FR')
%             text(0.2,0.9*ymaxval,sprintf('slope = %5.2f', c_FR(1) ),'Color','r','FontSize',14)
            title(sprintf('%s %s clu %i:  %i Hz FR response by period repetition',subject,session,clu,irate))
            
            % Save figure
            print_eps_kp(gcf,fullfile(savedir,sprintf('exUn_%s_%s_%i_clu%i_%iHz_FR',subject,session,channel,clu,irate)))
            
        end
        
        %--------- Delta Nspk ---------
        if iUn==ex_Un || p_wsr<0
            
            ymaxval = max(Nspk_paired(:))+1;
            
            figure; hold on
            plot([0 ymaxval],[0 ymaxval],'--k')
            plot(Nspk_paired(:,1),Nspk_paired(:,2),'.k','MarkerSize',30)
            if p_wsr<0.05
                plot(median(Nspk_paired(:,1)),median(Nspk_paired(:,2)),'.b','MarkerSize',50)
            end
            
            axis square
            xlim([0 ymaxval])
            ylim([0 ymaxval])
            if irate>2
                xlabel('nspks in first 250 ms')
                ylabel('nspks in last 250 ms')
            else 
                xlabel('nspks in first 500 ms')
                ylabel('nspks in last 500 ms')
            end
            text(0.02*ymaxval,0.95*ymaxval,sprintf('wsr p = %5.2f', p_wsr ),'Color','b','FontSize',14)
            title(sprintf('%s %s clu %i:  %i Hz N spks early vs late within trials',subject,session,clu,irate))
            
            % Save figure
            print_eps_kp(gcf,fullfile(savedir,sprintf('exUn_%s_%s_%i_clu%i_%iHz_dNspk',subject,session,channel,clu,irate)))
            
        end
        
        %------------ VS ------------
        if iUn==ex_Un || abs(c_VS(1))>100
            
            figure;
            plot(unique(x),VS,'.k','MarkerSize',40)
%             plotSpread(VS,'distributionIdx',unique(x),'distributionColors','k')
            hold on
            plot(unique(x),polyval(c_VS,unique(x)),'r-','LineWidth',2)
            
            ymaxval = 1;
            xlim([0 irep+1])
            ylim([0 ymaxval])
            xlabel('Repetition')
            ylabel('Vector Strength')
            text(0.2,0.9*ymaxval,sprintf('slope = %5.2f', c_VS(1) ),'Color','r','FontSize',14)
            title(sprintf('%s %s %i:  %i Hz Vector Strength by period repetition',subject,session,clu,irate))
            
            % Save figure
            print_eps_kp(gcf,fullfile(savedir,sprintf('exUn_%s_%s_%i_clu%i_%iHz_VS',subject,session,channel,clu,irate)))
            
        end
        
        %---------- Mean phase ----------
        if iUn==ex_Un || abs(c_Ph(1))>100
            
            figure;
            plot(unique(x),Ph,'.k','MarkerSize',40)
%             plotSpread(Ph,'distributionIdx',unique(x),'distributionColors','k')
            hold on
            plot(unique(x),polyval(c_Ph,unique(x)),'r-','LineWidth',2)
            
            ymaxval = 2*pi;
            xlim([0 irep+1])
            ylim([0 ymaxval])
            xlabel('Repetition')
            ylabel('Mean Phase')
            text(0.2,0.9*ymaxval,sprintf('slope = %5.2f', c_Ph(1) ),'Color','r','FontSize',14)
            title(sprintf('%s %s %i:  %i Hz Mean Phase by period repetition',subject,session,clu,irate))
            
            % Save figure
            print_eps_kp(gcf,fullfile(savedir,sprintf('exUn_%s_%s_%i_clu%i_%iHz_Phase',subject,session,channel,clu,irate)))
            
        end
        
    end %irate
    
    UnitData(iUn).SlopesFR      = Slopes_FR;
    UnitData(iUn).DeltaNspk     = DeltaNspk;
    UnitData(iUn).SlopesVS      = Slopes_VS;
    UnitData(iUn).SlopesPh      = Slopes_Ph;
    
    
end %iUn

AllSlopesFR  = vertcat(UnitData.SlopesFR);
AllDeltaNspk = vertcat(UnitData.DeltaNspk);
AllSlopesVS  = vertcat(UnitData.SlopesVS);
AllSlopesPh  = vertcat(UnitData.SlopesPh);


%% FR

hf_FR=figure;
plot([0 6],[0 0],'c','LineWidth',2)
hold on
plotSpread(AllSlopesFR,'distributionColors','k','showMM',3)
ylim([-2 2])
set(gca,'xtick',1:5,'xticklabel',Info.stim_ID_key(2:6))
ylabel('Slope')
title('Distribution of slopes: change in FR across repeated periods')

% Save figure
print_eps_kp(hf_FR,fullfile(savedir,'Slopes_FR'))


%% Delta N spikes

for irate = AMrates
    
    DeltaNspk_irate = vertcat(AllDeltaNspk{:,irate==AMrates});
    
    isig = DeltaNspk_irate(:,3)<0.05;
    ins  = DeltaNspk_irate(:,3)>=0.05;
    
    hf_dNS=figure;
    set(hf_dNS,'Position',tallfig)
    plot(DeltaNspk_irate(ins,1:2)','-k','LineWidth',1)
    hold on
    plot(DeltaNspk_irate(isig,1:2)','-b','LineWidth',5)
    xlim([0.75 2.25])
    ylim([0 15])%max(DeltaNspk_irate(:))+1])
    if irate>2
        set(gca,'xtick',1:2,'xticklabels',{'0-250 ms' '750-1000 ms'})
    else
        set(gca,'xtick',1:2,'xticklabels',{'0-500 ms' '500-1000 ms'})
    end
    ylabel('Number of spikes')
    title(sprintf('Change in n spikes early to late in %i Hz periodic',irate))
    
    % Save figure
    print_eps_kp(hf_dNS,fullfile(savedir,['DeltaNspks_' num2str(irate) 'Hz']))
    
end



%% VS

hf_VS=figure;
plot([0 6],[0 0],'c','LineWidth',2)
hold on
plotSpread(AllSlopesVS,'distributionColors','k','showMM',3)
ylim([-0.2 0.2])
set(gca,'xtick',1:5,'xticklabel',Info.stim_ID_key(2:6))
ylabel('Slope')
title('Distribution of slopes: change in VS across repeated periods')

% Save figure
print_eps_kp(hf_VS,fullfile(savedir,'Slopes_VS'))


%% Mean Phase

hf_Ph=figure;
plot([0 6],[0 0],'c','LineWidth',2)
hold on
plotSpread(AllSlopesPh,'distributionColors','k','showMM',3)
ylim([-1 1])
set(gca,'xtick',1:5,'xticklabel',Info.stim_ID_key(2:6))
ylabel('Slope')
title(sprintf('Distribution of slopes: change in resp Phase across repeated periods'))

% Save figure
print_eps_kp(hf_Ph,fullfile(savedir,'Slopes_Ph'))





end





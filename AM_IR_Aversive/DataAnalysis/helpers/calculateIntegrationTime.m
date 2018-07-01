function [IntegrationTime_spk, IntegrationTime_gap] = calculateIntegrationTime(spiketimes,TrialData,all_TDidx,Stream_Spikes,Stream_FRsmooth,AMrates,subject,session,channel,clu,RespType)

global fn

VSdata_spk    = nan(3,8);
VSdata_gap    = nan(3,8);
MeanPhase_spk = nan(1,8);
MeanPhase_gap = nan(1,8);
FR_raw        = nan(1,8);
GR_raw        = nan(1,8);

for istim = 2:6
    
    % Skip ITI stimuli
    ITIflag = 0;
    for iiti = 1:numel(ITIflag)
        
        period = 1000/AMrates(istim-1);
        
        TDidx = all_TDidx( TrialData.trID(all_TDidx)==istim & TrialData.ITIflag(all_TDidx) == ITIflag(iiti) );
        
        % Get timestamps of onsets and offsets
        clear t2 t3 Duration t_win
        t2 = TrialData.onset(TDidx);
        t3 = TrialData.offset(TDidx);
        Duration = mode(diff([t2 t3],1,2));
        t3 = t2 + Duration;
        
        % Preallocate
        raster = nan(numel(t2),Duration);
        PSTH   = nan(numel(t2),Duration);
        Spktimes = [];
        
        % Collect data for each trial
        nt=0; VStr_limit = 20;
        for it = randperm(numel(t2))
            nt=nt+1;
            
            % Collect spiking data (raw and normalized)
            raster(it,:) = Stream_Spikes( (t2(it)+1): t3(it) );
            PSTH(it,:) = Stream_FRsmooth( (t2(it)+1): t3(it) );
            
            % And collect spike data for VS calculation
            if istim>1 && istim<7  %&& nt<=VStr_limit
                sp=[]; sp = spiketimes( spiketimes>t2(it) & spiketimes<=t3(it) ) - t2(it);
                Spktimes = [Spktimes sp];
                % Check that 2 spike tracking methods match
                if ~all(find(raster(it,:))==sp)
                    keyboard
                end
            end
        end %it
        
        if any(any(isnan(raster)))
            keyboard
        end
        
        % Get mean FR for each stim
        FR_raw(1,istim) = mean(mean(raster,2)*1000);
        
        % Get VS data
        if (istim>1 && istim<7)
            
            % Calculate VS and mean phase for spikes
            [VSdata_spk(1,istim),VSdata_spk(2,istim),VSdata_spk(3,istim)] = vectorstrength(Spktimes,period);
            MeanPhase_spk(1,istim) = meanphase(sort(Spktimes),period);
            if MeanPhase_spk(1,istim)<0, keyboard, end
            
            %                         % Simulate Spktimes, to check if VS values match
            %                         Spktimes_sim_bin = poissrnd(repmat(mean(PSTH,1),size(PSTH,1),1)/1000);
            %                         Spktimes_sim=[];
            %                         for it=1:size(PSTH,1)
            %                             Spktimes_sim = [Spktimes_sim find(Spktimes_sim_bin(it,:))];
            %                         end
            %                         [VSdata_sim(1,istim),VSdata_sim(2,istim),VSdata_sim(3,istim)] = vectorstrength(Spktimes_sim,period);
            %                         MeanPhase_sim(1,istim) = meanphase(sort(Spktimes_sim),period);
            %
            
            % Simulate "Gaptimes", the inverted firing pattern
            PSTH_inv = -1*mean(PSTH,1) + min(mean(PSTH,1)) + max(mean(PSTH,1));
            Gaptimes_sim_bin = poissrnd(repmat(PSTH_inv,size(PSTH,1),1)/1000);
            Gaptimes_sim=[];
            for it=1:size(PSTH,1)
                Gaptimes_sim = [Gaptimes_sim find(Gaptimes_sim_bin(it,:))];
            end
            
            GR_raw(1,istim) = mean(mean(Gaptimes_sim_bin,2,'omitnan')*1000);
            
            % Calculate VS and mean phase for gaps (inverted
            % but balanced firing pattern)
            [VSdata_gap(1,istim),VSdata_gap(2,istim),VSdata_gap(3,istim)] = vectorstrength(Gaptimes_sim,period);
            MeanPhase_gap(1,istim) = meanphase(sort(Gaptimes_sim),period);
            if MeanPhase_gap(1,istim)<0, keyboard, end
        end
        
    end %iti
end %istim


%% Calculate integration time

% For SPIKES

% Wrap latency around to next cycle where necessary
PdcPhase_spk = MeanPhase_spk(2:6);
resets = 1+find(diff(PdcPhase_spk)<0);
for ii = 1:numel(resets)
    PdcPhase_spk(resets(ii):end) = PdcPhase_spk(resets(ii):end) + 2*pi;
end

% Fit slope to find integration time
weights_spk = FR_raw(1,2:6) .* VSdata_spk(1,2:6);

[IntegrationTime_spk,~,IT_mse] = lscov(AMrates',rad2deg(PdcPhase_spk'),weights_spk);

if any(rad2deg(PdcPhase_spk)<0)
    keyboard
end


% For GAPS

% Wrap latency around to next cycle where necessary
PdcPhase_gap = MeanPhase_gap(2:6);
resets = 1+find(diff(PdcPhase_gap)<0);
for ii = 1:numel(resets)
    PdcPhase_gap(resets(ii):end) = PdcPhase_gap(resets(ii):end) + 2*pi;
end

% Fit slope to find integration time
weights_gap = GR_raw(1,2:6) .* VSdata_gap(1,2:6);

[IntegrationTime_gap,~,IT_mse] = lscov(AMrates',rad2deg(PdcPhase_gap'),weights_gap);

if any(rad2deg(PdcPhase_gap)<0)
    keyboard
end


%% Plot response phases for this unit

hf_tmp=figure; hold on
ip(1)=plot(AMrates,rad2deg(PdcPhase_spk),'-k','LineWidth',2.5);
ip(2)=plot(AMrates,rad2deg(PdcPhase_gap),'-r','LineWidth',2.5);
for ir = 1:5
    plot(AMrates(ir),rad2deg(PdcPhase_spk(ir)),'.k','MarkerSize',10*weights_spk(ir))
    plot(AMrates(ir),rad2deg(PdcPhase_gap(ir)),'.r','MarkerSize',10*weights_gap(ir))
end
ylim([0 max( rad2deg(max(PdcPhase_spk)+pi/2), rad2deg(max(PdcPhase_gap)+pi/2)) ])
xlim([0 34])
axis square
%                 title(sprintf('Integration Time \nspks: %0.1f ms, \\color{red}gaps: %0.1f ms',IntegrationTime_spk,IntegrationTime_gap)...
%                     ,'Interpreter','tex')
title(sprintf('%s %s ch%i clu%i',subject,session,channel,clu),'Interpreter','none')
legend(ip,{sprintf('spks: %0.1f ms',IntegrationTime_spk) sprintf('gaps: %0.1f ms',IntegrationTime_gap)},'Location','southeast')
set(gca,'xtick',AMrates)
xlabel('AM rate (Hz)')
ylabel('Unwrapped mean phase (deg)')
hold off


% Save fig
savedir = fullfile(fn.processed,'Units',RespType);
if ~exist(savedir,'dir')
    mkdir(fullfile(savedir,'eps'))
    mkdir(fullfile(savedir,'svg'))
end
% print_eps_kp(hf_tmp,fullfile(savedir,'eps',sprintf('RespPhase_%s_%s_%i_%i',subject,session,channel,clu)))
% print_svg_kp(hf_tmp,fullfile(savedir,'svg',sprintf('RespPhase_%s_%s_%i_%i',subject,session,channel,clu)))





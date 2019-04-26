function IntegrationTime_spk = calculateIntegrationTime(spiketimes,TrialData,all_TDidx,Stream_Spikes,Stream_FRsmooth,AMrates,subject,session,channel,clu,RespType)

VSdata_spk    = nan(3,8);
MeanPhase_spk = nan(1,8);
FR_raw        = nan(1,8);

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
            
        end
        
    end %iti
end %istim


%% Calculate integration time

% To be valid, must have significant phase locking to 3 or more stimuli and
% a mean FR during pdc stim of at least 1 Hz!

if sum(VSdata_spk(3,2:6)<0.05)>2 && nanmean(FR_raw)>=1

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

fprintf('    integration time = %0.2f\n',IntegrationTime_spk)


%% Plot response phases for this unit

% hf_tmp=figure; hold on
% ip(1)=plot(AMrates,rad2deg(PdcPhase_spk),'-k','LineWidth',3);
% for ir = 1:5
%     plot(AMrates(ir),rad2deg(PdcPhase_spk(ir)),'.k','MarkerSize',15*weights_spk(ir)+1)
% end
% ylim([0 max( rad2deg(max(PdcPhase_spk)+pi/2) ) ])
% xlim([0 34])
% axis square
% title(sprintf('%s %s ch%i clu%i',subject,session,channel,clu),'Interpreter','none')
% legend(ip,sprintf('spks: %0.1f ms',IntegrationTime_spk),'Location','southeast')
% set(gca,'xtick',AMrates,'xscale','linear')
% xlabel('AM rate (Hz)')
% ylabel('Unwrapped mean phase (deg)')
% hold off


% % Save fig
% global fn
% savedir = fullfile(fn.processed,'Units',RespType);
% if ~exist(savedir,'dir')
%     mkdir(fullfile(savedir,'eps'))
%     mkdir(fullfile(savedir,'svg'))
% end
% print_eps_kp(hf_tmp,fullfile(savedir,'eps',sprintf('RespPhase_%s_%s_%i_%i',subject,session,channel,clu)))
% print_svg_kp(hf_tmp,fullfile(savedir,'svg',sprintf('RespPhase_%s_%s_%i_%i',subject,session,channel,clu)))


% close(hf_tmp)

else %not enough stim with sig VS 
    fprintf('    less than 3 periodic stim with significant VS\n')
    IntegrationTime_spk = nan;
%     IntegrationTime_gap = nan;
end


end


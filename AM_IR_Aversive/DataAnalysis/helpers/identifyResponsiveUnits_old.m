function [sigUnits,Resp] = identifyResponsiveUnits_old(Resp)
% sig_indices = identifyResponsiveUnits(Resp)
%
%  Identifies the units in the Resp input struct that are responsive to the
%  AM stimuli, according to any (or some combo) of:
%   - FR/trial distributions (periodic only)
%   - VS significant synchronization (periodic only)
%   - Rcorr classifier performance
%   - FR dprime formula (all stimuli compared to Silence)
% 
%  Next, add which stimuli are significant, and BMF for rate and VS.
%
%  KP, 2018-04  
% 


% Set thresholds
p_bonferonni_VS  =  0.0002;
% RC_M             =  3;

sigUnits = [];
for ii = 1:numel(Resp)
    
    % Rate tuning
    KW   =  Resp(ii).kw_p < 0.01;
    WX   =  Resp(ii).wx_p < 0.01;
    
    % Synchronization
    VS   =  any(Resp(ii).VSdata(3,:) < p_bonferonni_VS);
    
%     % Rcorr performance
%     RC_chance = 1/sum(~isnan(Resp(ii).FR_nrm));
%     RC   =  Resp(ii).RCorr.PC  >  (RC_M * RC_chance * 100);
    
    % IF EITHER FR OR VS IS SIGNIFICANT, KEEP UNIT
    if  (KW && WX)   ||   VS
        sigUnits = [sigUnits ii];
    end
    
    if  KW && WX
        % BMF - Rate
        FRpdc_tr = Resp(ii).FR_raw_tr(:,2:6);
        [~,Resp(ii).iBMF_FR] = max(mean(FRpdc_tr(sum(isnan(FRpdc_tr),2)==0,:),1));
    end
    
    if VS
        % BMF - Synchronization 
        BVS = max(Resp(ii).VSdata(1,(Resp(ii).VSdata(3,:) < p_bonferonni_VS)));
        Resp(ii).iBMF_VS = find(Resp(ii).VSdata(1,:)==BVS) - 1;
        % Significantly synchronized stimuli
        Resp(ii).iSync = find(Resp(ii).VSdata(3,:) < p_bonferonni_VS) - 1;
    end
    
end


end


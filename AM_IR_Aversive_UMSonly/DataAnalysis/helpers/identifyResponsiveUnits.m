function [sigUnits,UnitData] = identifyResponsiveUnits(UnitData)
% sig_indices = identifyResponsiveUnits(Resp)
%
%  Identifies the units in the Resp input struct that are responsive to the
%  AM stimuli, according to any (or some combo) of:
%   - FR/trial distributions (periodic only)
%   - VS significant synchronization (periodic only)
%   - Rcorr classifier performance
%   - FR dprime formula (all stimuli compared to Silence)
% 
%  Called by plotUnitDistributions, and re-saved.
%
%  KP, 2018-04 
% 


% Set thresholds
p_bonferonni_VS  =  0.0002;
% RC_M             =  3;

sigUnits = [];
for iUn = 1:numel(UnitData)
    
    % Rate tuning
    KW   =  UnitData(iUn).kw_p < 0.01;
    WX   =  UnitData(iUn).wx_p < 0.01;
    
    % Synchronization
    VS   =  any(UnitData(iUn).VSdata_spk(3,:) < p_bonferonni_VS);
    
%     % Rcorr performance
%     RC_chance = 1/sum(~isnan(Resp(ii).FR_nrm));
%     RC   =  Resp(ii).RCorr.PC  >  (RC_M * RC_chance * 100);
    
    % IF EITHER FR OR VS IS SIGNIFICANT, KEEP UNIT
    if  (KW && WX)   ||   VS
        sigUnits = [sigUnits iUn];
    end
    
    if  KW && WX
        % BMF - Rate
        FRpdc_tr = UnitData(iUn).FR_raw_tr(:,2:6);
        [~,UnitData(iUn).iBMF_FR] = max(mean(FRpdc_tr(sum(isnan(FRpdc_tr),2)==0,:),1));
    end
    
    if VS
        % BMF - Synchronization 
        BVS = max(UnitData(iUn).VSdata_spk(1,(UnitData(iUn).VSdata_spk(3,:) < p_bonferonni_VS)));
        UnitData(iUn).iBMF_VS = find(UnitData(iUn).VSdata_spk(1,:)==BVS) - 1;
        % Significantly synchronized stimuli
        UnitData(iUn).iSync = find(UnitData(iUn).VSdata_spk(3,:) < p_bonferonni_VS) - 1;
    end
    
end


end


function [AllCorrVals,TheseTrials] = checkCorrAct(CTTS,BSN,nSt,nUn,PickTrials,PSTHSize)

h=waitbar(0);

for iBS = 1:BSN
    
    waitbar(iBS/BSN,h,'Bootstrapping sets of training and testing trials...')
    
    
    %% If want all cells to have simultaneous trials, define all sets here
    
    if strcmp(PickTrials,'sim')
        %======================
        % SIMULTANEOUS trials
        %======================
        nTrials   = nan(nSt,1);
        
        PSTH_trs = nan(nSt,PSTHSize);
        
        for itmp = 1:nSt
            if ~isempty(find(~isnan(mean(CTTS(1,:,:,itmp),2)),1,'last'))
                nTrials(itmp) = find(~isnan(mean(CTTS(1,:,:,itmp),2)),1,'last');
            else
                keyboard
            end
            
            % TRAINING
            PSTH_trs(itmp,:)  = randperm(nTrials(itmp),PSTHSize);
            
        end
    end
    
    
    %% Preallocate
    
    AllCorrVals = cell(nSt,1);
    TheseTrials = cell(nSt,1);
    
    %% ###############################################################
    %                 Get rcorr values for each stim
    %  ###############################################################
    
    for ist = 1:nSt   %true stimulus of test trials
        
        SetRest = 1:nTrials(ist);
        SetRest(PSTH_trs(ist,:)) = [];
        
        CorrVals = nan(nUn,numel(SetRest));
        
        for iUn = 1:nUn
            
            T = mean(CTTS(iUn,:,PSTH_trs(ist,:),ist),3);
            
            for itr = 1:numel(SetRest)
                
                this_tr = SetRest(itr);
                
                % Get trial data
                trial_data = CTTS(iUn,:,this_tr,ist);
                
                % . .  . . . . . . . . . . . . . . . . . . . . . . . . . .
                % . Dot product of this test trial with each template
                if sum(trial_data)>0 && sum(T)>0
                    CorrVals(iUn,itr) = (trial_data * T') / (norm(trial_data)*norm(T));
                else
                    CorrVals(iUn,itr) = 0;
                end
            end %itr
        end %iUn
        
        AllCorrVals{ist} = CorrVals;
        TheseTrials{ist} = SetRest;
        
    end %ist
    
end %iBS

end
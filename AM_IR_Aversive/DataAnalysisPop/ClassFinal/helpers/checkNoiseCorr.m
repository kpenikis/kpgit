function checkNoiseCorr(Cell_Time_Trial_Stim)

figure;
xlim([-1 1])
ylim([-1 1])
hold on
xlabel('Signal correlation')
ylabel('Trial fidelity correlation')

% Find nTrials
nTrialMat = nan(size(Cell_Time_Trial_Stim,1),size(Cell_Time_Trial_Stim,4));
for ist = 1:size(Cell_Time_Trial_Stim,4)
    CT  = permute(sum(Cell_Time_Trial_Stim(:,:,:,ist),2),[1 3 2]);
    nTrialMat(:,ist) = sum(~isnan(CT),2);
end
Ntrs = min(nTrialMat,[],1);

Dur = size(Cell_Time_Trial_Stim,2);

% Calculate

SigCorrs = [];
TotCorrs = [];

for iC1 = 1:(size(Cell_Time_Trial_Stim,1)-1)
    for iC2 = (iC1+1):size(Cell_Time_Trial_Stim,1)
        
        % Gather data for this pair
        SignalC1 = [];
        SignalC2 = [];
        TotalC1  = [];
        TotalC2  = [];
        
        for iStim = 1:size(Cell_Time_Trial_Stim,4)
            
            SignalC1 = [SignalC1 mean(Cell_Time_Trial_Stim(iC1,:,1:Ntrs(iStim),iStim),3)];
            SignalC2 = [SignalC2 mean(Cell_Time_Trial_Stim(iC2,:,1:Ntrs(iStim),iStim),3)];
            
            % Raw signal on each trial
%             TotalC1 = [TotalC1 reshape(Cell_Time_Trial_Stim(iC1,:,1:Ntrs(iStim),iStim),1,Dur*Ntrs(iStim))];
%             TotalC2 = [TotalC2 reshape(Cell_Time_Trial_Stim(iC2,:,1:Ntrs(iStim),iStim),1,Dur*Ntrs(iStim))];
            
            % Trial FR variability
%             TotalC1 = [TotalC1 permute(mean(Cell_Time_Trial_Stim(iC1,:,1:Ntrs(iStim),iStim),2),[1 3 2])];
%             TotalC2 = [TotalC2 permute(mean(Cell_Time_Trial_Stim(iC2,:,1:Ntrs(iStim),iStim),2),[1 3 2])];
            
            % Fidelity of each trial to average 
            % dot product of stimulus mean activity with this trial
            Rcorrs1 = nan(1,Ntrs(iStim));
            Rcorrs2 = nan(1,Ntrs(iStim));
            for it = 1:Ntrs(iStim)  
                Rcorrs1(it) = mean(Cell_Time_Trial_Stim(iC1,:,1:Ntrs(iStim),iStim),3) * Cell_Time_Trial_Stim(iC1,:,it,iStim)';
                Rcorrs2(it) = mean(Cell_Time_Trial_Stim(iC2,:,1:Ntrs(iStim),iStim),3) * Cell_Time_Trial_Stim(iC2,:,it,iStim)';
            end
            TotalC1 = [TotalC1 Rcorrs1];
            TotalC2 = [TotalC2 Rcorrs2];
            
        end
        
        % Calculate signal correlation
        corrSig = corrcoef(SignalC1,SignalC2);
        SigCorrs =[SigCorrs corrSig(1,2)];
        
%         if corrSig(1,2)<-0.5
%             keyboard
%         end
        
        % Calculate total correlation
%         corrTot = corrcoef(TotalC1,TotalC2);
%         TotCorrs = [TotCorrs corrTot(1,2)];
        
        % Calculate total correlation
        corrTot = corrcoef(TotalC1,TotalC2);
        TotCorrs = [TotCorrs corrTot(1,2)];
        
    end
end


% figure;
plot(SigCorrs,TotCorrs,'.','MarkerSize',15)
% xlim([-1 1])
% ylim([-1 1])



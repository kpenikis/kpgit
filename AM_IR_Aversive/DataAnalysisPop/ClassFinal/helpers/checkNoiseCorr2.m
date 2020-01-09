function checkNoiseCorr2(Cell_Time_Trial_Stim,CR)

% Find nTrials
nTrialMat = nan(size(Cell_Time_Trial_Stim,1),size(Cell_Time_Trial_Stim,4));
for ist = 1:size(Cell_Time_Trial_Stim,4)
    CT  = permute(sum(Cell_Time_Trial_Stim(:,:,:,ist),2),[1 3 2]);
    nTrialMat(:,ist) = sum(~isnan(CT),2);
end
Ntrs = min(nTrialMat,[],1);

% Calculate
SigCorrs = [];
for iC1 = 1:(size(Cell_Time_Trial_Stim,1)-1)
    for iC2 = (iC1+1):size(Cell_Time_Trial_Stim,1)
        
        % Gather data for this pair
        SignalC1 = [];
        SignalC2 = [];
        
        for iStim = 1:size(Cell_Time_Trial_Stim,4)
            
            SignalC1 = [SignalC1 mean(Cell_Time_Trial_Stim(iC1,:,1:Ntrs(iStim),iStim),3)];
            SignalC2 = [SignalC2 mean(Cell_Time_Trial_Stim(iC2,:,1:Ntrs(iStim),iStim),3)];
                        
        end
        
        % Calculate signal correlation
        corrSig = corrcoef(SignalC1,SignalC2);
        SigCorrs =[SigCorrs corrSig(1,2)];
        
    end
end

pcdiffs = CR.PC(2:2:end)-CR.PC(1:2:end);
dpdiffs = CR.dprime(2:2:end)-CR.dprime(1:2:end);


% DPRIME
figure;
plot(SigCorrs,dpdiffs,'.','MarkerSize',15)

% Fit line to data 
x = linspace(-1,1,100);
c = polyfit(SigCorrs',dpdiffs,1);
disp(['Equation is y = ' num2str(c(1)) '*x + ' num2str(c(2))])
y_est = polyval(c,x);
% Add trend line to plot
hold on
plot(x,y_est,'r--','LineWidth',2)

xlim([-1 1])
ylabel('Sim - Shuff d''')

[r1,p1]=corr(SigCorrs',dpdiffs)


% PC
figure;
scatter(SigCorrs,pcdiffs,30,CR.PC(2:2:end),'filled')

% Fit line to data 
x = linspace(-1,1,100);
c = polyfit(SigCorrs',pcdiffs,1);
disp(['Equation is y = ' num2str(c(1)) '*x + ' num2str(c(2))])
y_est = polyval(c,x);
% Add trend line to plot
hold on
plot(x,y_est,'r--','LineWidth',2)

xlim([-1 1])
ylabel('Sim - Shuff PC')

[r2,p2]=corr(SigCorrs',pcdiffs)




% Load data
load('IntTimeData_kp')

% Preallocate
IT_result = nan(size(RespPhase_deg,1),1);
IT_mse    = nan(size(RespPhase_deg,1),1);
figure; hold on

% for each unit
for iu = 1:size(RespPhase_deg,1)
    
    % Plot response phase vs AM rate 
    plot(AMrates,RespPhase_deg(iu,:),'-')
    
    % Calculate integration time
    [IT_result(iu,1), ~, IT_mse(iu,1)] = lscov( AMrates', RespPhase_deg(iu,:)', Ntrials(iu,:) );
    % could weight with several different things, including:   
    %   Ntrials(iu,:)     
    %   VS(iu,:)   
    %   VS_pval(iu,:).^-1
    
end

% View distribution of Integration Times across population
figure;
histogram(IT_result,50)




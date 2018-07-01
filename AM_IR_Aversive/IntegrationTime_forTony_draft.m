

AMrates       = [2 4 8 16 32];
RespPhase_deg = vertcat(Resp.Phase_deg);
Ntrials       = vertcat(Resp.IT_ntr);

VS      = nan(size(Ntrials));
VS_pval = nan(size(Ntrials));

for iu = 1:numel(Resp)
    VS(iu,:)       = vertcat(Resp(iu).VSdata(1,2:6));
    VS_pval(iu,:)  = vertcat(Resp(iu).VSdata(3,2:6));
end

fn = set_paths_directories('','',1);
save(fullfile(fn.processed,'IntTimeData_kp'),'AMrates','RespPhase_deg','Ntrials','VS','VS_pval','-v7.3')



clear all

fn = set_paths_directories('','',1);

% Load data
load(fullfile(fn.processed,'IntTimeData_kp'))

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




function runALL_anPlots(subject,session)
% runALL_anPlots(subject,session)
%   Plot responses of all clusters of a session to all stimuli, for all
%   response metrics.
%   
%   Inputs
%     subject: subject name as string
%     session: session label as string
%     channel: (optional) channel number as double
%         clu: (optional) cluster label as double
% 
% KP, 02-2017,
% 



% Load Data structure
datadir  = '/Users/kpenikis/Documents/SanesLab/Data/AMJitter/ProcessedData';
filename = sprintf( '%s_sess-%s_Data',subject,session);
load(fullfile(datadir,subject,filename));


METRICS = {'FR' 'FF' 'FF-avPds' 'FF-Pds' 'VS' 'VS-Pds' 'RS' 'RS-Pds' 'standardFR' 'Corr'};

% Go through each cluster and call analysis programs
for channel = 1:numel(Data.ch)
    for iu = 1:numel(Data.ch(channel).clu)
        
        clu = Data.ch(channel).clu(iu).label(1);
        
        close all
                
        for im = 1:numel(METRICS)
            ap_xJitter(subject,session,channel,clu,METRICS{im})
        end
        
        
    end 
end 

end
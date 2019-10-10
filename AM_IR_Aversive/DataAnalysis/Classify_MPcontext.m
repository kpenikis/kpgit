function Data = Classify_MPcontext(USE_MEASURE,SHUFFLE,ii)
%
%  Classify_MPcontext(USE_MEASURE,SHUFFLE_SPIKES,ii)
% 
%   After splitting up ClassifyMPHcontext_matched, this program runs
%   classifier on matched MP data struct that already exists.
%
%  KP, 2019-09
%


global fn AMrates Iterations WinLen

%!!!!!!!!!!!!!!!!!!!
if nargin<1
    USE_MEASURE = 'FR'; 'spikes';
end
%!!!!!!!!!!!!!!!!!!!
Iterations =  10000;
%!!!!!!!!!!!!!!!!!!!
SHUFFLE    =  0;
%!!!!!!!!!!!!!!!!!!!
RERUN      =  1;
%!!!!!!!!!!!!!!!!!!!
WinLen     =  62;
%!!!!!!!!!!!!!!!!!!!
theseRates = 1:4;
%!!!!!!!!!!!!!!!!!!!


%% Load data files

fn = set_paths_directories;

q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q
%-------
spkshift = mean([UnitData([UnitData.IntTime_spk]>0).IntTime_spk]);
%-------

% Load data from matched MPs
q = load(fullfile(fn.processed,'MatchedMPs','MatchedMPdata'));
Data = q.Data;
clear q


if ~isstring(WinLen)
    WinLen_str = num2str(WinLen);
end


% Assign savename according to shuffled or not
savedir = fullfile(fn.processed,'MPcontext');
savename = ['MPcontextSU_' WinLen_str ];acti

% Preallocate output
if RERUN
    Res = struct;
    Un1 = 1;
else
    load(fullfile(savedir,'tmp',savename));
    Un1 = size(Res,1)+1;
end


%% Run classifier 

for iUn = Un1:size(Data,1)
    
    fprintf('Un %i running classifier... ',iUn)
    
    if SHUFFLE
        keyboard
        % shuffle trials/timebins across stimuli/contexts
    end
    
    for irate = theseRates 
        
        if ~isfield(Data(iUn,irate).data,'Context')
            continue
        end
        
        % Run classifier and save results
        fprintf('%iHz ',AMrates(irate))
        Res(iUn,irate).L1o  = get_classifier_data( Data(iUn,irate).data, -1 );
        
    end  
    
    %% Save Res to tmp folder
    
    save(fullfile(savedir,'tmp',savename),'Res','-v7.3')
    fprintf('\n')
    
    
end %iUn

% Save final Data struct
save(fullfile(savedir,savename),'Res','-v7.3')


end %function





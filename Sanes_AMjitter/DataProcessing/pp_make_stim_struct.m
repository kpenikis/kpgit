
function [Stim, epoc_trs] = pp_make_stim_struct( epData, block, stimfiles )
%
%  pp_make_stim_struct( epData, block )
%    Sub-routine called by pp_prepare_format.
%    
%    Creates matlab struct with number of elements equal to number of
%    trials, and fields with stimulus information for each trial. 
%    
%    This function should be modified for each experimenter/experiment, to
%    include the relevant parameters and source of trial onset timestamps.
%    
%  KP, 2016-04; last updated 2016-10
% 
 

% STIM ONSET TIMES & PARAMS


if isfield(epData.epocs,'rvID')  %AM jitter recording
    fileID  = epData.epocs.rvID.data';
    AMrate  = epData.scalars.Pars.data(2,:);
    AMdepth = epData.epocs.Dpth.data';
    AMonset = round( 1000*  epData.scalars.Pars.ts );

    onsets  = round( 1000*  epData.epocs.rvID.onset )';  %ms
    stimDur = round( 1000* (epData.epocs.rvID.offset - epData.epocs.rvID.onset) )';  %ms
    dB_level = epData.epocs.dBSP.data';
    
    HP      = epData.scalars.Pars.data(3,:);
    LP      = epData.scalars.Pars.data(4,:);
    
    Freq    = nan(size(onsets));
    
elseif isfield(epData.epocs,'AMrt') %AM and level tuning
    
    AMrate  = epData.epocs.AMrt.data';
    AMdepth = epData.epocs.Dpth.data';
    AMonset = round( 1000*  epData.scalars.Pars.ts );

    onsets  = round( 1000*  epData.epocs.AMrt.onset )';  %ms
    stimDur = 10+round( 1000* (epData.epocs.AMrt.offset - epData.epocs.AMrt.onset) )';
    dB_level = epData.epocs.dBSP.data';
    
    HP      = epData.scalars.Pars.data(2,:);
    LP      = epData.scalars.Pars.data(3,:);
    if all(HP==21) && all(LP==21)
        HP = 100*ones(size(HP));
        LP = 20000*ones(size(HP));
    end
    
    Freq    = nan(size(onsets));
    fileID  = nan(size(onsets));
    
elseif isfield(epData.epocs,'Freq') %Frequency tuning
    Freq    = epData.epocs.Freq.data;
    
    onsets  = round( 1000*  epData.epocs.Freq.onset )';  %ms
    stimDur = 5+round( 1000* (epData.epocs.Freq.offset - epData.epocs.Freq.onset) )';
    dB_level = epData.scalars.Pars.data(1,:);
        
    HP      = nan(size(onsets));
    LP      = nan(size(onsets));
    fileID  = nan(size(onsets));
    AMrate  = nan(size(onsets));
    AMdepth = nan(size(onsets));
    AMonset = nan(size(onsets));

end

%% Mark trials as 0: passive uncontrolled, 1: behaving, 2: passive drinking.
% Code will always stop here for jitter stim blocks, and require user input
% about subject's drinking behavior.
N = 0; 
if all(~strcmp(stimfiles,'basic_tuning'))
    N = get_drinking_trials(block);
end
if N==0
    behaving = zeros(numel(onsets),1);
elseif isnan(N)
    behaving = 2*ones(1,numel(onsets));
else
    behaving = [2*ones(1,N) zeros(1,numel(onsets)-N)];
end
%%

% Check that the 2 saved sources from openex correspond
if length(AMonset) ~= length(onsets)
    keyboard
    %%%check trial params saved in both ways
% %     AMonset - onsets; stimDur;
% %     epoc_trs = [1:44 46:48];
% %     %%%keep only full trials
% %     %epoc pars
% %     onsets   = onsets(epoc_trs);
% %     dB_level = dB_level(epoc_trs);
% %     stimDur  = stimDur(epoc_trs);
% %     Freq     = Freq(epoc_trs);
% %     AMdepth  = AMdepth(epoc_trs);
% %     fileID   = fileID(epoc_trs);
    
    %par pars
else
    epoc_trs = 1:length(onsets);
end


% LOAD BEHAVIORAL DATA



% CREATE STIM STRUCT

Stim = struct();
for it = 1:numel(onsets)
    
    % Basic info from epData struct
    Stim(it).onset       =    onsets(it);
    Stim(it).onsetAM     =    AMonset(it);
    Stim(it).behaving    =    behaving(it);
    Stim(it).block       =    block;
    
    Stim(it).dB          =    dB_level(it);
    Stim(it).stimDur     =    stimDur(it);
    
    Stim(it).Freq        =    Freq(it);
    
    Stim(it).HP          =    HP(it);
    Stim(it).LP          =    LP(it);
    
    Stim(it).AMrate      =    AMrate(it);
    Stim(it).AMdepth     =    AMdepth(it);
    
    Stim(it).fileID      =    fileID(it);
    
    if iscell(stimfiles)
        Stim(it).stimfile    =    stimfiles{fileID(it)};
    else
        Stim(it).stimfile    =    stimfiles;
    end
    
    
    % incorporate behavioral responses here later
    
end


end

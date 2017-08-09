
function [Stim, epoc_trs] = pp_make_stim_struct( epData, block )
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


if isfield(epData.epocs,'AMrt') %AM and level tuning
    
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
    
elseif isfield(epData.epocs,'Freq') %Frequency tuning
    Freq    = epData.epocs.Freq.data;
    
    onsets  = round( 1000*  epData.epocs.Freq.onset )';  %ms
    stimDur = 5+round( 1000* (epData.epocs.Freq.offset - epData.epocs.Freq.onset) )';
    dB_level = epData.scalars.Pars.data(1,:);
        
    HP      = nan(size(onsets));
    LP      = nan(size(onsets));

end


%%

epoc_trs = 1:length(onsets);


% LOAD BEHAVIORAL DATA



% CREATE STIM STRUCT

Stim = struct();
for it = 1:numel(onsets)
    
    % Basic info from epData struct
    Stim(it).onset       =    onsets(it);
    Stim(it).block       =    block;
    
    Stim(it).dB          =    dB_level(it);
    Stim(it).stimDur     =    stimDur(it);
    
    Stim(it).Freq        =    Freq(it);
    
    Stim(it).HP          =    HP(it);
    Stim(it).LP          =    LP(it);
    
        
    % incorporate behavioral responses here 
    
end


end

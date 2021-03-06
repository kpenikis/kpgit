function [subject, session, channel, Spikes] = pp_launch_manual_sort( subject, session, channel, Spikes )
%
%  pp_launch_manual_sort( subject, session, channel )
%    Basically just launches the splitmerge tool from UMS. But if called
%    with output arugments, will add all relevant variable to workspace,
%    which facilitates updating and saving of the Spikes file at the end of
%    manual sorting. ie. the output args here are the input args of
%    pp_save_manual_sort function. 
% 
%  KP, 2016-04; last updated 2017-04
% 

% Set directories based on processor used 
fn = set_paths_directories;

clear spikes

% Load data structures

fprintf('loading data...\n')
if ~exist('Spikes','var')
    filename = sprintf('%s_sess-%s_Spikes',subject,session);
    load(fullfile(fn.processed,subject,filename));
end
filename = sprintf('%s_sess-%s_Info',subject,session);
load(fullfile(fn.processed,subject,filename));


% If a manual sort has already been saved, open from where left off
% Also serves as a quick way to start over - set ch flag back to 0. 

if Spikes.man_sort(channel)==1 
    Spikes.sorted(channel).params.display.default_waveformmode = 2;
    
    splitmerge_tool(  Spikes.sorted(channel)  )
    
elseif Spikes.man_sort(channel)==0 
    Spikes.clustered(channel).params.display.default_waveformmode = 2;
    
    splitmerge_tool(  Spikes.clustered(channel)  )
    
end



% Add automatic screening for good clusters
% clus = Spikes.clustered(channel).labels(:,1);
% for ic = 1:numel(clus)
%     [ev,lb,ub,RPV] = ss_rpv_contamination( Spikes.clustered(channel), clus(ic)  )
% end
% Output:
%   ev   - expected value of % contamination,
%   lb   - lower bound on % contamination, using 95% confidence interval
%   ub   - upper bound on % contamination, using 95% confidence interval
%   RPV  - number of refractory period violations observed in this cluster



end
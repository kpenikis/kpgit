function cludata = format_add_Data( ib, ip, id, blocks, bk_vec,...
                            LP_HP_dB_rate, depths, jitter, behav,...
                            fieldname, datastruct, cludata, Waves )
% cludata = format_add_Data( ib, ip, id, blocks, bk_vec, LP_HP_dB_rate, depths, stimfn, shift, r_mean, p_mean, r_trials )
%   Formats basic analysis data to be added to the Data structure.
%   No use in checking if stim values already filled, because all analyses
%   called from the same code, within the same loops, with the same checks.
% 

cludata.stim(blocks==ib).pars(ip).block                 = bk_vec;
cludata.stim(blocks==ib).pars(ip).pars                  = LP_HP_dB_rate(ip,:);
cludata.stim(blocks==ib).pars(ip).depths(id).depth      = depths(id);
cludata.stim(blocks==ib).pars(ip).depths(id).jitter     = jitter;
cludata.stim(blocks==ib).pars(ip).depths(id).behav      = behav;
if ~isempty(Waves)
    cludata.stim(blocks==ib).pars(ip).depths(id).Waves  = Waves;
end

cludata.stim(blocks==ib).pars(ip).depths(id).(fieldname) = datastruct;


end
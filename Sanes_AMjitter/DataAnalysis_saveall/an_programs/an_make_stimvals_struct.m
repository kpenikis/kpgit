function stim = an_make_stimvals_struct( id, depths, jitter, behav, stim )
% cludata = format_add_Data( ib, ip, id, blocks, bk_vec, LP_HP_dB_rate, depths, stimfn, shift, r_mean, p_mean, r_trials )
%   Formats basic analysis data to be added to the Data structure.
%   No use in checking if stim values already filled, because all analyses
%   called from the same code, within the same loops, with the same checks.
% 

stim(id).depth      = depths(id);
stim(id).jitter     = jitter;
stim(id).behav      = behav;

end
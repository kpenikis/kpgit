function add_datapoints_table(stimvals, dp_struct, session, cluname, clulabel, stimpars, behav, indVar)

global DataTable

stimvals = stimvals(strcmp(behav,stimvals(:,3)),:);

for istim = 1:size(stimvals,1)
    
    this_jitter = char(stimvals(istim,1));
    this_depth  = char(stimvals(istim,2));
    this_behav = strcmp(behav,{dp_struct.classifFR0.behav});
    
    switch indVar
        case 'depth'
            cond_idx = find(strcmp(this_jitter,dp_struct.stim(:,1)));
            ind_val  = convert_depth_proptodB(str2num(this_depth));
            whichind = 1;
        case 'jitter'
            cond_idx = find(strcmp([this_depth '_depth'],dp_struct.stim(:,1)));
            ind_val  = str2num(strtok(this_jitter,'_'));
            whichind = find(strcmp(this_jitter, stimvals( strcmp(num2str(ind_val), strtok(stimvals(:,1),'_')) & strcmp(this_depth,stimvals(:,2)) ,1) ));
    end
    try
    if (strcmp(indVar,'jitter') && ind_val==0) || (strcmp(indVar,'depth') && ind_val==convert_depth_proptodB(0))
        dp_FR_cls = nan;
        dp_FR_fml = nan;
        dp_SpV100 = nan;
        dp_SpV25  = nan;
        dp_SpV10  = nan;
        dp_Corr   = nan;
    else
        dp_FR_cls = dp_struct.classifFR0(this_behav).dprime.data{cond_idx}     (round(dp_struct.classifFR0(this_behav).dprime.data{cond_idx}(:,1) - ind_val,4)==0 ,2);
        %%%%%
        dp_FR_fml = dp_struct.formulaFR.output(this_behav).dprime_mat{cond_idx}(round(dp_struct.formulaFR.output(this_behav).dprime_mat{cond_idx}(:,1) - ind_val,4)==0 ,2);
        %%%%%
        dp_SpV100 = dp_struct.classifSpV100(this_behav).dprime.data{cond_idx}  (round(dp_struct.classifSpV100(this_behav).dprime.data{cond_idx}(:,1) - ind_val,4)==0 ,2);
        dp_SpV25  = dp_struct.classifSpV25(this_behav).dprime.data{cond_idx}   (round(dp_struct.classifSpV25(this_behav).dprime.data{cond_idx}(:,1)  - ind_val,4)==0 ,2);
        dp_SpV10  = dp_struct.classifSpV10(this_behav).dprime.data{cond_idx}   (round(dp_struct.classifSpV10(this_behav).dprime.data{cond_idx}(:,1)  - ind_val,4)==0 ,2);
        dp_Corr   = dp_struct.classifCorr0(this_behav).dprime.data{cond_idx}   (round(dp_struct.classifCorr0(this_behav).dprime.data{cond_idx}(:,1)  - ind_val,4)==0 ,2);

    end
    
    %'Session' 'cluname' 'SU/MU' 'HP LP' 'dBSPL' 'CenterRate' 'BehavState' 'Depth' 'Jitter' 'dprime'... 'indVar'
    DT_addrow = {session cluname ...
        discretize(clulabel, [0.5 1.5 2.5 3.5 4.5], 'categorical',{'unk', 'SU', 'MU', 'noise'})...
        [stimpars(1) stimpars(2)] stimpars(3) stimpars(4) ...
        categorical(behav, {'P' 'D' 'A'}, {'Passive' 'Drinking' 'Active Behavior'})...
        this_depth this_jitter dp_FR_cls(whichind) dp_SpV100(whichind) dp_SpV25(whichind) dp_SpV10(whichind) dp_Corr(whichind) dp_FR_fml(whichind) indVar};
    
    DataTable = [DataTable; DT_addrow];
    catch
        keyboard
    end
    
end

end




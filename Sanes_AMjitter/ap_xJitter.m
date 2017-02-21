function ap_xJitter(subject,session,channel,clu,METRIC)
% ap_xJitter(subject,session,channel,clu,METRIC)
%   Calculates the response of this cluster to each stimulus, reflected by
%   the response measure specified by input variable METRIC. Plots and 
%   saves as figures, and inserts each datapoint into the Data file. 
% 
%   This program combines ap_barplot_indJitter, ap_heatplot_indJitter, and
%   ap_timeseries_indJitter. 
% 
%   If runing through all clusters of a session, this program is called by 
%   the program runALL_anPlots.
% 
%   Inputs
%     subject: subject name as string
%     session: session label as string
%     channel: channel number as double
%         clu: cluster label as double
%      METRIC: string that serves as code for desired analyses measure
%              possible values for METRIC are:
%              {'FR' 'FF' 'FF-avPds' 'FF-Pds' 'VS' 'VS-Pds' 'RS' 'RS-Pds' 'standardFR' 'Corr'}
% 
% KP, 02-2017.
% 


set(0,'DefaultAxesFontSize',10)
set(0,'DefaultTextInterpreter','none')

addpath(genpath('analysis_modules'),'helpers','plotting_modules')


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Load raster struct
[raster,Data] = get_raster(subject,session,channel,clu);

% Remove stimuli with fewer than 8 trials
raster = raster(cellfun(@length,{raster.tr_idx}) > 8);


% Set up cludata struct
labels = vertcat(Data.ch(channel).clu.label);
iclu = labels(:,1)==clu;

if isfield(Data.ch(channel).clu(iclu),'analyses')
    cludata = Data.ch(channel).clu(iclu).analyses;
else 
    cludata = struct;
end


% Get blocks for loop and designate which ones to combine
[blocks,group_blocks] = set_grouped_blocks(unique([raster.block]));

% Make separate figures for each param except jitter and depth
for ib = blocks
    
    % Get data for these blocks
    bk_raster = raster([raster.block]==ib);
    
    % Add data from blocks set to combine, if needed
    if ~isempty(group_blocks(:,1)==ib)
        bk_raster = [bk_raster raster([raster.block]== group_blocks(group_blocks(:,1)==ib,2) )];
        bk_str = [num2str(ib) num2str(group_blocks(group_blocks(:,1)==ib,2))];
        bk_vec = [ib group_blocks(group_blocks(:,1)==ib,2)];
    else
        bk_str = num2str(ib);
        bk_vec = ib;
    end
    
    % Find unique stimuli based on other parameters
    [LP_HP_dB_rate,~,np] = unique([bk_raster.HP; bk_raster.LP; bk_raster.dB; bk_raster.AMrate]','rows');
    
    
    for ip = 1:max(np)
        
        param_raster = bk_raster(np==ip);
        
        [depths,~,ndpth] = unique([param_raster.AMdepth]);
        
        % Calculate baseline rate for this condition
        baselineFR = calc_baselineFR(param_raster);
        
        hF = figure; 
        
        % Go through by increasing AMdepth
        for id = 1:max(ndpth)
            
            % Collapse data across blocks, if params repeated
            stim = collapse_blocks(param_raster(ndpth==id));
            
            
            %% GET DATAPOINT
            % 
            Waves=[];
            try
            switch METRIC
                
                case 'Corr'
                    [r_mean,p_mean,shift,r_trials,Waves] = corr_spks(stim,subject);
                    
                    % Format output to add to Data struct
                    datastruct.shift = shift;
                    datastruct.rSt   = r_mean;
                    datastruct.pSt   = p_mean;
                    datastruct.rTr   = r_trials;
                    fieldname = 'Corr';
                    
                case 'psth'
                    keyboard
%                     [r_mean,p_mean,shift] = corr_psth(stim,binsize,subject);

                case 'FR'
                    [data_mean,data_std,data] = calc_FR(stim);
                    
                    % Format output to add to Data struct
                    datastruct.mean  = data_mean';
                    datastruct.std   = data_std;
                    datastruct.trs   = data;
                    fieldname = 'FR';
                    
                case 'FF'
                    binsize = 250;
                    [data_mean,data_std] = calc_FF(stim,binsize);
                    data = nan;
                    
                    % Format output to add to Data struct
                    datastruct.mean    = data_mean;
                    datastruct.prctl   = data_std;
                    datastruct.binsize = binsize;
                    fieldname = 'FF';
                    
                case {'FF-avPds' 'FF-Pds'}
                    [data] = calc_FF_periods(stim,subject);
                    data_mean = mean(data,2,'omitnan');
                    data_std  = [std(data,0,2,'omitnan') std(data,0,2,'omitnan')]./2;
                    
                    % Format output to add to Data struct
                    datastruct.mean  = data_mean;
                    datastruct.stdPd = data_std;
                    datastruct.trs   = data;
                    fieldname = 'FFpds';
                    
                case {'VS' 'VS-Pds'}
                    [data_mean,data] = calc_VSRS(stim,subject);
                    data_std = nan(length(data_mean),2);
                    
                    % Format output to add to Data struct
                    datastruct.mean    = data_mean;
                    datastruct.periods = data;
                    fieldname = 'VS';
                    
                case {'RS' 'RS-Pds'}
                    [~,~,data_mean,data] = calc_VSRS(stim,subject);
                    data_std = nan(length(data_mean),2); 
                    
                    % Format output to add to Data struct
                    datastruct.mean    = data_mean;
                    datastruct.periods = data;
                    fieldname = 'RS';
                    
                case 'standardFR'
                    [data_mean,data_std,data] = calc_standardFR(stim,subject);
                    
                    % Format output to add to Data struct
                    datastruct.mean  = data_mean';
                    datastruct.std   = data_std;
                    datastruct.trs   = data;
                    fieldname = 'standardFR';

            end
            catch
                keyboard
            end
            
            
            %% Add datapoint to cludata
            
            cludata = format_add_Data( ib, ip, id, blocks, bk_vec,...
                                        LP_HP_dB_rate, depths, [stim.jitter], {stim.behaving},...
                                        fieldname, datastruct, cludata, Waves );
            
            
            %% PLOT DATAPOINT
            
            figure(hF)
            nsp = [ (1+sum(ndpth<id)) : (sum(ndpth<id)+sum(ndpth==id)) ];
            switch METRIC
                case {'Corr' 'FR' 'FF' 'FF-avPds' 'VS' 'RS' 'standardFR'}
                    subplot(1,length(ndpth), nsp ,'align')
                case {'FF-Pds' 'VS-Pds' 'RS-Pds'}
                    subplot(length(ndpth), 1, nsp ,'align')
            end
            hold on
            pause(0.1)

            % Call subplot function
            switch METRIC
                
                case 'Corr'
                    subplot_ts(stim,r_mean,p_mean,shift)
                    pause(0.4)
                    
                case {'FR' 'FF' 'FF-avPds' 'VS' 'RS' 'standardFR'}
                    subplot_bar(stim,data_mean,data_std,data,baselineFR,METRIC);
                    
                case {'FF-Pds' 'VS-Pds' 'RS-Pds'}
                    lL = max(nsp)==length(ndpth);
                    subplot_mtx(stim,data,METRIC,lL);
                    
            end
            
            
        end
        
        
        %% FORMAT PLOT
        
        % Set axis limits to be same for all subplots
        hAllAxes = findobj(hF,'type','axes');
        
        switch METRIC
            
            case 'Corr'
                if iscell(get(hAllAxes,'YLim'))
                    ylims = [min(cellfun(@min,get(hAllAxes,'YLim'))) max(cellfun(@max,get(hAllAxes,'YLim')))];
                else
                    ylims = get(hAllAxes(1),'YLim');
                end
                set(hAllAxes,'YLim',ylims)
                
                
            case {'FR' 'FF' 'FF-avPds' 'VS' 'RS' 'standardFR'}
                
                hAllAxes = findobj(hF,'type','axes');
                if iscell(get(hAllAxes,'YLim'))
                    ymax = max(cellfun(@max,get(hAllAxes,'YLim')));
                else
                    ymax = max(get(hAllAxes(1),'YLim'));
                end
                set(hAllAxes,'YLim',[0 ymax])
                
        end
        
        % Get string for savename
        str_pars = strsplit(num2str(LP_HP_dB_rate(ip,:)));
        if numel(depths)==1
            str_dpth = num2str(depths*100);
        else
            str_dpth = 'VARIED';
        end
        
        % Add stimulus info in title
        title_str = sprintf('ch %i clu%i\nnoise: %s - %s Hz  |  %sdB\nAM %s Hz  |  %s dpth  |  blk%s',...
            channel,clu,str_pars{1},str_pars{2},str_pars{3},str_pars{4},str_dpth,bk_str);
        suptitle(title_str);
        
        
        
        %% SAVE FIGURE
        
        % Specify correct subfolder and savename
        switch METRIC
            case 'Corr'
                savefolder = 'Corr_sp-stim';
                savename   = sprintf('%s_%s_CorrSpkAM_ch%i_clu%i_AM%sHz_%sdpth_%sdB_%s-%s_blk%s',...
                    subject,session,channel,clu,str_pars{4},str_dpth,str_pars{3},str_pars{1},str_pars{2},bk_str);
            case 'FR'
                savefolder = 'FRavg';
                savename   = sprintf('%s_%s_FR_ch%i_clu%i_AM%sHz_%sdpth_%sdB_%s-%s_blk%s',...
                    subject,session,channel,clu,str_pars{4},str_dpth,str_pars{3},str_pars{1},str_pars{2},bk_str);
            case 'FF'
                savefolder = ['FFavg-' num2str(binsize)];
                savename   = sprintf('%s_%s_FFbin%i_ch%i_clu%i_AM%sHz_%sdpth_%sdB_%s-%s_blk%s',...
                    subject,session,binsize,channel,clu,str_pars{4},str_dpth,str_pars{3},str_pars{1},str_pars{2},bk_str);
            case 'FF-avPds'
                savefolder = 'FFavg-periods';
                savename   = sprintf('%s_%s_FFpdAVG_ch%i_clu%i_AM%sHz_%sdpth_%sdB_%s-%s_blk%s',...
                    subject,session,channel,clu,str_pars{4},str_dpth,str_pars{3},str_pars{1},str_pars{2},bk_str);
            case 'FF-Pds'
                savefolder = 'FF-periods';
                savename   = sprintf('%s_%s_FFpd_ch%i_clu%i_AM%sHz_%sdpth_%sdB_%s-%s_blk%s',...
                    subject,session,channel,clu,str_pars{4},str_dpth,str_pars{3},str_pars{1},str_pars{2},bk_str);
            case 'VS'
                savefolder = 'VSavg';
                savename   = sprintf('%s_%s_VS_ch%i_clu%i_AM%sHz_%sdpth_%sdB_%s-%s_blk%s',...
                    subject,session,channel,clu,str_pars{4},str_dpth,str_pars{3},str_pars{1},str_pars{2},bk_str);
            case 'VS-Pds'
                savefolder = 'VS-periods';
                savename   = sprintf('%s_%s_VSpd_ch%i_clu%i_AM%sHz_%sdpth_%sdB_%s-%s_blk%s',...
                    subject,session,channel,clu,str_pars{4},str_dpth,str_pars{3},str_pars{1},str_pars{2},bk_str);
            case 'RS'
                savefolder = 'RSavg';
                savename   = sprintf('%s_%s_RS_ch%i_clu%i_AM%sHz_%sdpth_%sdB_%s-%s_blk%s',...
                    subject,session,channel,clu,str_pars{4},str_dpth,str_pars{3},str_pars{1},str_pars{2},bk_str);
            case 'RS-Pds'
                savefolder = 'RS-periods';
                savename   = sprintf('%s_%s_RSpd_ch%i_clu%i_AM%sHz_%sdpth_%sdB_%s-%s_blk%s',...
                    subject,session,channel,clu,str_pars{4},str_dpth,str_pars{3},str_pars{1},str_pars{2},bk_str);
            case 'standardFR'
                savefolder = 'standardFR';
                savename   = sprintf('%s_%s_standardFR_ch%i_clu%i_AM%sHz_%sdpth_%sdB_%s-%s_blk%s',...
                    subject,session,channel,clu,str_pars{4},str_dpth,str_pars{3},str_pars{1},str_pars{2},bk_str);
                
        end
        

        
        
        % Set save directory 
        datadir  = '/Users/kpenikis/Documents/SanesLab/Data/AMJitter/ProcessedData';
        an_dir = fullfile(datadir,subject,'^an_plots',session,savefolder);
        if ~exist(an_dir,'dir')
            mkdir(an_dir)
        end
        
        % Save figure
        set(gcf,'PaperOrientation','landscape');
        print(hF,'-dpdf',fullfile(an_dir,savename),'-bestfit')
%         print(hF,'-depsc',fullfile(an_dir,savename))
        
        
        
    end
end


%% Add to Data struct and save file

labels = vertcat(Data.ch(channel).clu.label);
iclu = labels(:,1)==clu;
Data.ch(channel).clu(iclu).analyses = cludata;
filename = sprintf( '%s_sess-%s_Data',subject,session); 
save(fullfile(datadir,subject,filename),'Data','-v7.3');



end


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



function ap_timeseries_indJitter(subject,session,channel,clu,METRIC)
% Calculates various response measures for each stimulus, as defined in input
% variable METRIC, and plots the results as bar plots.


set(0,'DefaultAxesFontSize',10)
set(0,'DefaultTextInterpreter','none')

addpath(genpath('analysis_modules'),'helpers','plotting_modules')

%~~~~~~~~~~~~~~~~~~~~
% Load raster struct
[raster,Data] = get_raster(subject,session,channel,clu);

% labels = vertcat(Data.ch(channel).clu.label);
% iclu = labels(:,1)==clu;
% if isfield(Data.ch(channel).clu(iclu),'analyses')
%     cludata = Data.ch(channel).clu(iclu).analyses;
% else 
    cludata = struct;
% end


% Remove stimuli with fewer than 8 trials
raster = raster(cellfun(@length,{raster.tr_idx}) > 8);


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
        
        hFw = figure; %timeseries
%         hFt = figure; %barplot 
        hold on
        
        % Go through by increasing AMdepth
        for id = 1:max(ndpth)
            
            % Collapse data across blocks, if params repeated
            stim = collapse_blocks(param_raster(ndpth==id));
            
            % Get data to plot
            try
            switch METRIC
                case 'spks'
                    [r_mean,p_mean,shift,r_trials] = corr_spks(stim,subject);
                case 'psth'
                    keyboard
%                     [r_mean,p_mean,shift] = corr_psth(stim,binsize,subject);
            end
            catch
                keyboard
            end
            
            
            % Format output to add to Data struct
            corr.shift = shift;
            corr.rSt   = r_mean;
            corr.pSt   = p_mean;
            corr.rTr   = r_trials;

            cludata = format_add_Data( ib, ip, id, blocks, bk_vec,...
                                        LP_HP_dB_rate, depths, {stim.stimfn},...
                                        'corr', corr, cludata );
            
            
%             % Plot subplot timeseries
%             figure(hFw)
%             nsp = [ (1+sum(ndpth<id)) : (sum(ndpth<id)+sum(ndpth==id)) ];
%             subplot(1,length(ndpth), nsp ,'align')
%             hold on, 
%             pause(0.1)
%             % Call subplot function
%             subplot_ts(stim,r_mean,p_mean,shift)
%             pause(0.4)
            
            % Plot subplot barplot
% %             figure(hFt)
% %             nsp = [ (1+sum(ndpth<id)) : (sum(ndpth<id)+sum(ndpth==id)) ];
% %             subplot(1,length(ndpth), nsp ,'align')
% %             hold on
% %             pause(0.1)
% %             % Call subplot function
% %             subplot_bar(stim,r_trials,nan(length(r_trials),2),nan,[],METRIC);
% %             pause(0.4)


        end
        
        
        
        
        % Set axis limits to be same for all subplots
        hAllAxes = findobj(hFw,'type','axes');
        if iscell(get(hAllAxes,'YLim'))
            ylims = [min(cellfun(@min,get(hAllAxes,'YLim'))) max(cellfun(@max,get(hAllAxes,'YLim')))];
        else
            ylims = get(hAllAxes(1),'YLim');
        end
        set(hAllAxes,'YLim',ylims)
        
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
        
        
        % Save figure
        
        switch METRIC
            case 'spks'
                savefolder = 'Corr_sp-stim';
                savename   = sprintf('%s_%s_FRresp_ch%i_clu%i_AM%sHz_%sdpth_%sdB_%s-%s_blk%s',...
                    subject,session,channel,clu,str_pars{4},str_dpth,str_pars{3},str_pars{1},str_pars{2},bk_str);
            case 'psth'
                savefolder = 'Corr_psth-stim';
                savename   = sprintf('%s_%s_FF-bin%i_ch%i_clu%i_AM%sHz_%sdpth_%sdB_%s-%s_blk%s',...
                    subject,session,binsize,channel,clu,str_pars{4},str_dpth,str_pars{3},str_pars{1},str_pars{2},bk_str);
        end
        
        datadir  = '/Users/kpenikis/Documents/SanesLab/Data/AMJitter/ProcessedData';
        an_dir = fullfile(datadir,subject,'^an_plots',session,savefolder);
        if ~exist(an_dir,'dir')
            mkdir(an_dir)
        end
        
        set(gcf,'PaperOrientation','landscape');
        print(hFw,'-dpdf',fullfile(an_dir,savename),'-bestfit')
%         print(hF,'-depsc',fullfile(an_dir,savename))
        
        
        
    end
end

keyboard

% Add to Data sturct and save file
Data.ch(channel).clu(Data.ch(channel).clu.label(1)==clu).analyses = cludata;
filename = sprintf( '%s_sess-%s_Data',subject,session); 
save(fullfile(datadir,subject,filename),'Data','-v7.3');

end


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



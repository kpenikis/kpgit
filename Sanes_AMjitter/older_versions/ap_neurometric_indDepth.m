function ap_neurometric_indDepth(subject,session,channel,clu,METRIC,merge_blocks)
% ap_neurometric_indDepth
%   runs neurometric classifier on spike data


%%%%%%%%%%%%%%%%%%%%%%%%
% IMPORTANT PARAMETER! %
binsize = 50;
%%%%%%%%%%%%%%%%%%%%%%%%


set(0,'DefaultAxesFontSize',10)
set(0,'DefaultTextInterpreter','none')

% IF SAVING PDF FILES
figFontSize      = 16;
setMarkerSize    = 30;
setLineWidth     = 3;

% IF SAVING EPS FILES
% figFontSize      = 24;
% rasterMarkerSize = 18;
% rasterLineWidth  = 1;

addpath(genpath('analysis_modules'),'helpers')


% Load raster struct
raster = get_raster(subject,session,channel,clu);

% Remove stimuli with fewer than 8 trials
raster = raster(cellfun(@length,{raster.tr_idx}) > 8);


if ~merge_blocks
    % Get blocks for loop and designate which ones to combine
    [blocks,group_blocks] = set_grouped_blocks(unique([raster.block]));
else
    % this case if you do want to merge across blocks (i.e. depth
    % discrimination when blocks of 75% depth and 41% depth are presented
    % separately). 
    keyboard
    blocks = 1;
end

% Make separate figures for each param except jitter and depth
for ibk = blocks
    
    if ~merge_blocks
        % Get data for these blocks
        bk_raster = raster([raster.block]==ibk);
        
        % Add data from blocks set to combine, if needed
        if ~isempty(group_blocks(:,1)==ibk)
            bk_raster = [bk_raster raster([raster.block]== group_blocks(group_blocks(:,1)==ibk,2) )];
            bk_str = [num2str(ibk) num2str(group_blocks(group_blocks(:,1)==ibk,2))];
        else
            bk_str = num2str(ibk);
        end
        
    else
        bk_raster = raster;
        bk_str = 'all';
    end


% Find unique stimuli based on other parameters
[LP_HP_dB_rate,~,np] = unique([bk_raster.HP; bk_raster.LP; bk_raster.dB; bk_raster.AMrate]','rows');


for ip = 1:max(np)
    
    param_raster = bk_raster(np==ip);
    
    % Separate figures by behavioral state
    for ib = {'D' 'P' 'A'}
        
        stim = param_raster(strcmp({param_raster.behaving},ib));
        
        % Get unique depths and jitter vectors
        [depths,~,ndpth] = unique([stim.AMdepth]);
        [fIDs,~,nfid] = unique({stim.stimfn},'stable');
        
        if numel(depths)<3  %skip data that can't yield a depth function
            continue
        end
        
        mdur = min([stim.stimDur]);
        dp = nan( numel(fIDs), numel(depths) );
        
        % Get nogo stimulus info
        NGidx = ndpth==find(depths==0);
        NOGO = stim(NGidx); 
        if numel(NOGO)>1
            NOGO = collapse_blocks(NOGO);
        elseif numel(NOGO)<1
            keyboard
        end
        
        % Get nogo spiketimes, starting when sound begins
        nogo_x = NOGO.x(NOGO.x>0);
        nogo_y = NOGO.y(NOGO.x>0);
        
        % Step through the go stimuli to compare to nogo
        for fid = 1:max(nfid)
            for id = 1:max(ndpth)
                if numel(stim((nfid==fid)&(ndpth==id)))<1, continue, end
                
                GO = collapse_blocks(stim((nfid==fid)&(ndpth==id)));
                
                % Get spiketimes for go stimulus, starting when sound begins
                try
                go_x = GO.x(GO.x>0);
                go_y = GO.y(GO.x>0);
                catch
                    keyboard
                end
                
                % Clip spiketimes at end of stimulus                
                nogo_y = nogo_y(nogo_x<=mdur);
                nogo_x = nogo_x(nogo_x<=mdur);
                go_y   = go_y(go_x<=mdur);
                go_x   = go_x(go_x<=mdur);
                
                % Set binsize depending on decoder type
                switch METRIC
                    case 'FR'
                        bin = mdur;
                    case 'spT'
                        bin = binsize;
                end
                
                % Convert raster format to accomodate distance calculations
                GOs = zeros(max(go_y),mdur);
                rsGO = zeros(max(go_y),floor(mdur/bin));
                for iy = 1:max(go_y)
                    GOs(iy,go_x(go_y==iy)) = 1;
                    rsGO(iy,:) = sum( reshape( GOs(iy,1:(bin*floor(mdur/bin))) , [bin,floor(mdur/bin)] ) ,1);
                end
                
                NOGOs = zeros(max(nogo_y),mdur);
                rsNOGO = zeros(max(nogo_y),floor(mdur/bin));
                for iy = 1:max(nogo_y)
                    NOGOs(iy,nogo_x(nogo_y==iy)) = 1;
                    rsNOGO(iy,:) = sum( reshape( NOGOs(iy,1:(bin*floor(mdur/bin))) , [bin,floor(mdur/bin)] ),1);
                end
                
                
                %%
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                % Run classifier
                dp(fid,id) = run_classifier(rsGO,rsNOGO,1000);
                
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                
                
            end
        end
        
        % If depth is 0 (unmodulated), use same d' value for both
        % jitter conditions, since it's irrelevant
        if sum(~isnan(dp(:,depths==0)))==1
            dp(isnan(dp(:,depths==0)),depths==0) = dp(~isnan(dp(:,depths==0)),depths==0);
        elseif any(isnan(dp(:,depths==0)))
            keyboard
        end
        
        % Convert depth % to dB re 100
        depths(depths==0) = 0.01; depths(depths==1) = 0.99;
        depth_dB = 20 .* log(depths(2:end));
        
        
        %~~~~~~~~~~~~~~~~~~~~~~
        % Fit data with weibull
% %         try
% %         bestFitParams=nan(max(nfid),2);
% %         for fid = 1:max(nfid)
% %             bestFitParams(fid,:) = fit_weibull(depth_dB,dp(fid,2:end));
% %         end
% %         catch
% %             keyboard
% %         end
        
        %~~~~~~~~~~~~~~~~~~~~~~
        %   Plot everything
        pcols = copper(max(nfid));
        
        hF=figure;
        hP=plot( repmat([20 .* log(depths(1)) depth_dB],max(nfid),1)' , dp' ,...
            '.','MarkerSize',setMarkerSize);
        for ih = 1:numel(hP)
            hP(ih).Color = pcols(ih,:);
            if any(strcmp(fIDs{ih},{'AM_4Hz_0.mat' 'AM_8Hz_0.mat' 'AM_16Hz_0.mat'}))
                hP(ih).Color = 'blue';
            end
        end
        hold on
        plot([-100 0],[1 1],'--k')
        xlim([-100 0])
        ylim([0 3.5])
        set(gca,'XTick',[-100 -10 -1])
        xlabel('AM depth (dB re 100%)')
        hL=legend(fIDs,'Location','northwest'); hL.Interpreter = 'none';
        set(gca,'FontSize',figFontSize)
        
% %         % Plot best fit weibull function over observed data range
% %         plot_stim = -42:0.25:0;
% %         for fid = 1:max(nfid)
% %             Weibull_fit = real(weibull(bestFitParams(fid,1),bestFitParams(fid,2),plot_stim-min(depth_dB),1));
% %             wp= plot(plot_stim,Weibull_fit,'-','LineWidth',setLineWidth, 'Color', pcols(fid,:));
% %             if any(strcmp(fIDs{fid},{'AM_4Hz_0.mat' 'AM_8Hz_0.mat' 'AM_16Hz_0.mat'}))
% %                 wp.Color = 'blue';
% %             end
% %         end
        
        hold off
        
        % Add stimulus info in title
        str_pars = strsplit(num2str(LP_HP_dB_rate(ip,:)));
        title_str = sprintf('ch %i clu%i\nnoise: %s - %s Hz  |  %sdB\nAM %s Hz',...
            channel,clu,str_pars{1},str_pars{2},str_pars{3},str_pars{4});
        
        switch METRIC
            case 'FR'
                ylabel('d-prime, based on FR');
                t=title([title_str '  --  FR']);
            case 'spT'
                ylabel('d-prime, based on spike vector');
                t=title([title_str '  --  spike vector, binsize ' num2str(bin)]);
        end
        
        % Set font sizes
        set(findall(gcf, 'Type','text'), 'FontSize', figFontSize)
        t.FontSize=0.75*figFontSize;
        hold off
        
        
        %~~~~~~~~~~~~~~
        % Save figure %
        %~~~~~~~~~~~~~~
        
        switch METRIC
            case 'FR'
                savefolder = 'Classifier-FR';
                savename   = sprintf('%s_%s_dPrime_FR_ch%i_clu%i_AM%sHz_%sdB_%s-%s_blk%s',...
                    subject,session,channel,clu,str_pars{4},str_pars{3},str_pars{1},str_pars{2},bk_str);
            case 'spT'
                savefolder = ['Classifier-spT' ];
                savename   = sprintf('%s_%s_dPrime_spT_bin%i_ch%i_clu%i_AM%sHz_%sdB_%s-%s_blk%s',...
                    subject,session,bin,channel,clu,str_pars{4},str_pars{3},str_pars{1},str_pars{2},bk_str);
        end
        
        datadir  = '/Users/kpenikis/Documents/SanesLab/Data/AMJitter/ProcessedData';
        an_dir = fullfile(datadir,subject,'^an_plots',session,savefolder);
        if ~exist(an_dir,'dir')
            mkdir(an_dir)
        end
        
        set(gcf,'PaperOrientation','landscape');
        print(hF,'-dpdf',fullfile(an_dir,savename),'-bestfit')
%         print(hF,'-depsc',fullfile(an_dir,savename))
        
    end
end

end











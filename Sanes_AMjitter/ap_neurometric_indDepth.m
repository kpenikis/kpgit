function ap_neurometric_indDepth(subject,session,channel,clu,METRIC,merge_blocks,raster)


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

addpath('analysis_modules','helpers')

% Load raster struct for this unit
if nargin<6 || ~exist('raster','var')
    savedir  = '/Users/kpenikis/Documents/SanesLab/Data/processed_data';
    savename = sprintf('%s_sess-%s_raster_ch%i_clu%i',subject,session,channel,clu);
    load(fullfile(savedir,subject,savename))
end


% Remove stimuli with fewer than 8 trials
raster = raster(cellfun(@length,{raster.tr_idx}) > 8);


if ~merge_blocks
    % Find blocks and designate which ones to combine
    blocks = unique([raster.block]);
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    group_blocks = [89 90];  %only 2 at a time for now
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for ic = size(group_blocks,1)
        blocks(blocks==group_blocks(ic,2)) = [];
    end
else
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
                
                
                % Run classifier
                dp(fid,id) = run_classifier(go_x,go_y,nogo_x,nogo_y,1000);
                if dp(fid,id)>2
                    keyboard
                end
                
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
        depth_dB = 20 .* log(depths);
        
        pcols = copper(max(nfid));
        
        hF=figure;
        hP=plot( repmat(depth_dB,fid,1)' , dp' ,...
            '.-','LineWidth',setLineWidth,'MarkerSize',setMarkerSize);
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
        ylabel('d-prime, based on FR')
        xlabel('AM depth (dB re 100%)')
        hL=legend(fIDs,'Location','northwest'); hL.Interpreter = 'none';
        set(gca,'FontSize',figFontSize)
        hold off
        
        % Add stimulus info in title
        str_pars = strsplit(num2str(LP_HP_dB_rate(ip,:)));
        title_str = sprintf('ch %i clu%i\nnoise: %s - %s Hz  |  %sdB\nAM %s Hz',...
            channel,clu,str_pars{1},str_pars{2},str_pars{3},str_pars{4});
        t=title(title_str);
        
        % Set font sizes
        set(findall(gcf, 'Type','text'), 'FontSize', figFontSize)
        t.FontSize=0.75*figFontSize;
        hold off
        
        
        %~~~~~~~~~~~~~~
        % Save figure %
        %~~~~~~~~~~~~~~
        an_dir = fullfile(savedir,subject,'^an_plots',session);
        if ~exist(an_dir,'dir')
            mkdir(an_dir)
        end
        
        switch METRIC
            case 'FR'
                savename = sprintf('%s_%s_dPrime_FR_ch%i_clu%i_AM%sHz_%sdB_%s-%s_blk%s',...
                    subject,session,channel,clu,str_pars{4},str_pars{3},str_pars{1},str_pars{2},bk_str);
            case 'spT'
                savename = sprintf('%s_%s_dPrime_spT_ch%i_clu%i_AM%sHz_%sdB_%s-%s_blk%s',...
                    subject,session,binsize,channel,clu,str_pars{4},str_dpth,str_pars{3},str_pars{1},str_pars{2},bk_str);
        end
        
        set(gcf,'PaperOrientation','landscape');
        print(hF,'-dpdf',fullfile(an_dir,savename),'-bestfit')
%         print(hF,'-depsc',fullfile(an_dir,savename))
        
    end
end

end











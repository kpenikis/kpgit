function dp = an_responses(METRIC,rasters,indVar,subject,plotOptions,ALLcolors,xlimits,PLOT_DIFF)
% an_responses(METRIC,rasters,indVar,parname,which_classifier,figHandle)
%   Calculates the response of this cluster to each stimulus, reflected by
%   the response measure specified by input variable METRIC. Plots and
%   saves as figures, and inserts each datapoint into the Data file.
%   Program is based on an_xJitter, but is more generalized.
%
%   If running for all units, called by updatePopulationPlots.
%
%   Inputs
%      METRIC: string that serves as code for desired analyses measure
%              possible values for METRIC are:
%              {'FR' 'FF' 'FF-avPds' 'FF-Pds' 'VS' 'VS-Pds' 'RS' 'RS-Pds' 'standardFR' 'Corr'}
%
% KP, 04-2017.
%
% global session unit

% Prepare variables for output
dp_depth  = [];
dp_jitter = [];
dp_data   = [];
dp_baseFR = [];

% Get independent and condition variable values
switch indVar
    case 'jitter'
        indVar_vals = unique([rasters.jitter]);
        condVar = 'AMdepth';
        condVar_vals = unique([rasters.AMdepth]);
        condVar_vals(condVar_vals==0) = []; %if jitter on x axis, leave out depth of 0
    case 'depth'
        indVar = 'AMdepth';
        indVar_vals = unique([rasters.AMdepth]);
        condVar = 'jitter';
        condVar_vals = unique([rasters.jitter]);
end

% Calculate baseline rate for this condition
baselineFR = calc_baselineFR(rasters);


% Set up empty matrix
data_plot = nan(numel(unique([rasters.(indVar)])),numel(condVar_vals));

% Go through each condition
for ic = 1:numel(condVar_vals)
    
    clear xvals
    
    switch indVar
        case 'jitter'
            cond_rasters = rasters([rasters.(condVar)]==condVar_vals(ic));
            icondVar{ic} = [num2str(condVar_vals(ic)) '_depth'];
            xvals = str2double(strtok([cond_rasters.(indVar)],'_'))';
            check_x = str2double(strtok(unique([rasters.(indVar)]),'_'))';
        case 'AMdepth'
            cond_rasters = rasters(strcmp([rasters.(condVar)],condVar_vals{ic}));
            icondVar{ic} = condVar_vals{ic};
            xvals = convert_depth_proptodB([cond_rasters.(indVar)])';
            check_x = convert_depth_proptodB(unique([rasters.(indVar)])');
    end
    
    
    %% GET DATAPOINTS
    %
    try
        switch METRIC
            case 'maxCorr'
                % Get correlations at all lags for all stimuli
                [Rs,~,sh,~,~,~] = corr_spks(cond_rasters,subject);
                
                % Go through each stimulus and find max corr at its best
                % lag value.
                for ii = 1:size(Rs,1)
                    
                    [~,peakLag]=findpeaks(Rs(ii,:),'MinPeakProminence',0.005);
                    if numel(peakLag)<1
                        [~,peakLag]=findpeaks(Rs(ii,:),'MinPeakProminence',0.003);
                    end
                    if numel(peakLag)<1
                        disp('weak correlations. setting lag to 0')
                        peakLag = find(sh==0);
                    elseif numel(peakLag)>1
                        [~,leastshift] = min(abs(peakLag-find(sh==0)));
                        peakLag = peakLag(leastshift);
                    end
                    
                    [~,~,~,Rtrs,~,~]   = corr_spks(cond_rasters(ii),subject,sh(peakLag));
                    
                    data(ii,1) = mean(Rtrs);
                    
                end
                
                data_plot(ismember(check_x,xvals,'rows'),ic) = data;
                
            case 'shftCorr'
                % Get correlations at all lags for all stimuli
                [Rs,~,sh,~,~,~] = corr_spks(cond_rasters,subject);
                
                % Go through each stimulus and find its best lag value.
                for ii = 1:size(Rs,1)
                    
                    [~,peakLag]=findpeaks(Rs(ii,:),'MinPeakProminence',0.005);
                    if numel(peakLag)<1
                        [~,peakLag]=findpeaks(Rs(ii,:),'MinPeakProminence',0.003);
                    end
                    if numel(peakLag)<1
                        disp('weak correlations. setting lag to 0')
                        peakLag = find(sh==0);
                    elseif numel(peakLag)>1
                        [~,leastshift] = min(abs(peakLag-find(sh==0)));
                        peakLag = peakLag(leastshift);
                    end
                                        
                    data(ii,1) = sh(peakLag);
                
                end
                
                data_plot(ismember(check_x,xvals,'rows'),ic) = data;
                
            case 'FR'
                [data,data_std] = calc_FR(cond_rasters);
                data_plot(ismember(check_x,xvals,'rows'),ic) = data;%(data-baselineFR)./baselineFR;
                
            case 'FF'
                binsize = 250;
                data = calc_FF(cond_rasters,binsize);
                data_plot(ismember(check_x,xvals,'rows'),ic) = data;
                
                if any(data<0)
                    keyboard
                end
                
            case {'FFavPds' 'FF_Pds'}
                [data_all] = calc_FF_periods(cond_rasters,subject);
                data = mean(data_all,2,'omitnan');
                data_plot(ismember(check_x,xvals,'rows'),ic) = data;
                data_std  = [std(data,0,2,'omitnan') std(data,0,2,'omitnan')]./2;
                
            case {'VS' 'VS_Pds'}
                data = calc_VSRS(cond_rasters,subject);
                data_plot(ismember(check_x,xvals,'rows'),ic) = data;
                data_std = nan(length(data),2);
                
                if any(data>1)
                    keyboard
                end
                
            case {'RS' 'RS_Pds'}
                [~,~,data,~] = calc_VSRS(cond_rasters,subject);
                data_plot(ismember(check_x,xvals,'rows'),ic) = data;
                data_std = nan(length(data),2);
                
            case 'standardFR'
                [data,data_std] = calc_standardFR(cond_rasters,subject);
                data_plot(ismember(check_x,xvals,'rows'),ic) = data;
        end
    catch
        keyboard
    end
    
    
    % Collect datapoints for table output
    dp_depth  = [dp_depth; [cond_rasters.AMdepth]' ];
    dp_jitter = [dp_jitter; [cond_rasters.jitter]'];
    dp_data   = [dp_data; data ];
    dp_baseFR = [dp_baseFR; repmat(baselineFR,size(data)) ];
    
    
    
end %ic

% Apploy same value to both jitters for AM depth of 0
if strcmp(indVar,'AMdepth') && isnan(data_plot(1,1))
    data_plot(1,1) = data_plot(1,2);
end

% Plot datapoint
try
    if ~PLOT_DIFF
        for ic = 1:numel(condVar_vals)
            plot(xvals,data_plot(:,ic),'.-','MarkerSize',30,...
                'Color', ALLcolors( strcmp(strtok(icondVar{ic},'_'),strtok(plotOptions.colSelect,'_')), : ), ...
                'LineWidth', 2)
        end
        ylabel(METRIC)
        
    elseif PLOT_DIFF && numel(condVar_vals)==2 && strcmp(icondVar{1},'0')
        plot(xvals,diff(data_plot,1,2),'.-','MarkerSize',30,...
            'Color', ALLcolors( strcmp(strtok(icondVar{2},'_'),strtok(plotOptions.colSelect,'_')), : ), ...
            'LineWidth', 2)
        plot([min(xlimits) max(xlimits)],[0 0],'--k','LineWidth',0.75)
        ylabel([METRIC ': jitter cond - periodic cond'])
        
    else
        keyboard
    end
catch
    keyboard
end

set(gca,'xlim',xlimits)


dp = table(dp_depth,dp_jitter,dp_data,dp_baseFR);
dp.Properties.VariableNames = {'depth' 'jitter' METRIC 'baselineFR' };

end


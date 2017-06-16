function fullRaster_periodic( rasters, session, cluname, clulabel,...
                        stimpars, behav, indVar )
                    

% global fn
% 
% binsize = 25; 
% 
% rasters = rasters([rasters.AMdepth]==0.75);
% if isempty(rasters), return, end
% 
% Jitters = [rasters.jitter];
% 
% baseFR = calc_baselineFR(rasters);
% 
% 
% % Set options for plot
% [~, plotOptions] = setOptions;
% plotOptions.colSelect = {'0_' '10_' '20_' '30_' '40_' '50_' '70_' '100_' '150_' '200_'};
% ALLcolors = copper( numel(plotOptions.colSelect) );
% set(0,'DefaultTextInterpreter','none')
% set(0,'DefaultAxesFontSize',16)
% 
% nRsp = ceil((numel(Jitters))/2);
% 
% scrsz = get(0,'ScreenSize');
% if nRsp>1
%     figsize = [1 scrsz(4) scrsz(3)/2 scrsz(4)];
% else
%     figsize = [1 scrsz(4)/2.5 scrsz(3)/2 scrsz(4)/2.5];
% end
% 
% % Create figures
% 
% hf = figure; 
% set(gcf,'Position',figsize,'NextPlot','add');
% axis tight
% hold on
% 
% 
% 
% %% First get periodic stimulus data
% 
% idxPdc = find(strcmp([rasters.jitter],'0'));
% if idxPdc~=1, keyboard, end
% 
% try
% % Prep for other periodic pd spikes
% rV = jitter_LUT(rasters(idxPdc).AMrate,char(rasters(idxPdc).jitter));
% tV = cumsum( [ ceil(rasters(idxPdc).AMonset+0.75*1000/rV(1)) 1000./rV(2:end) ]);
% catch
%     keyboard
% end
% 
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% clear x y
% [x,y,prev250FR,prev100FR] = standardPd_getspikes(rasters(idxPdc));
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
% x_p=[]; y_p=[]; p_p=[];
% these_pds = [2 6];
% for ipd = 1:numel(these_pds)
%     
%     this_pd = these_pds(ipd);
%     
%     x_p = [x_p rasters(idxPdc).x( (rasters(idxPdc).x >= tV(this_pd-1)) & (rasters(idxPdc).x < tV(this_pd)) ) - tV(this_pd-1) + 1 ];
%     y_p = [y_p rasters(idxPdc).y( (rasters(idxPdc).x >= tV(this_pd-1)) & (rasters(idxPdc).x < tV(this_pd)) ) ];
%     p_p = [p_p repmat(this_pd, size(rasters(idxPdc).x( (rasters(idxPdc).x >= tV(this_pd-1)) & (rasters(idxPdc).x < tV(this_pd)) )) )];
%     
% end
% 
% spikes  = zeros(max(y),250);
% spikes2 = zeros(max(y),250);
% spikes6 = zeros(max(y),250);
% 
% for iy = 1:max(y)
%     
%     spikes(iy,x(y==iy)) = 1;
%     
%     spikes2(iy, x_p(y_p==iy & p_p==2) ) = 1;
%     spikes6(iy, x_p(y_p==iy & p_p==6) ) = 1;
%     
% end
% 
% nTr_pdc = max(y);
% 
% 
% %......................................................................
% 
% % Calculate smoothed PSTH
% FRsmooth_STANDARD = smoothPSTH_v2(spikes);
% SDsmooth_STANDARD = smooth_STD(spikes,binsize);
% 
% % [FRsmooth_STANDARD,~,~,~,SDbin_STANDARD,SDbinT_STANDARD] = smoothPSTH_v2(spikes);
% FRsmooth_2 = smoothPSTH_v2(spikes2);
% SDsmooth_2 = smooth_STD(spikes2,binsize);
% FRsmooth_6 = smoothPSTH_v2(spikes6);
% SDsmooth_6 = smooth_STD(spikes6,binsize);
% 
% 
% 
% %% Now plot PSTH of periods 2 and 6 on top of standard
% 
% % figure(hf(1))
% % subplot(nRsp,2,1); hold on
% % title(['periodic, ntr=' num2str(numel(rasters(idxPdc).tr_idx))])
% % for it = 1:(max(y)+max(y_p(p_p==2))+max(y_p(p_p==6)))
% %     plot([0 250],[it it],'Color',[0.5 0.5 0.5],'LineWidth', 0.5*plotOptions.rasterLineWidth)
% % end
% % plot(x_p(p_p==2),max(y)+y_p(p_p==2),'r+','MarkerSize',plotOptions.markerSize, 'LineWidth', plotOptions.rasterLineWidth)
% % plot(x,y,'k+','MarkerSize',plotOptions.markerSize, 'LineWidth', plotOptions.rasterLineWidth)
% % plot(x_p(p_p==6),max(y)+max(y_p(p_p==2))+y_p(p_p==6),'b+','MarkerSize',plotOptions.markerSize, 'LineWidth', plotOptions.rasterLineWidth)
% % ylim([0 max(y)+max(y_p(p_p==2))+max(y_p(p_p==6))+1])
% % 
% % figure(hf(2))
% % subplot(nRsp,2,1); hold on
% % title(['periodic, ntr=' num2str(numel(rasters(idxPdc).tr_idx))])
% % plot(FRsmooth_STANDARD,'k','LineWidth',plotOptions.lineWidth)
% % plot(FRsmooth_2,'r','LineWidth',plotOptions.lineWidth)
% % plot(FRsmooth_6,'b','LineWidth',plotOptions.lineWidth)
% 
% figure(hf)
% subplot(nRsp,2,1); hold on
% title(['periodic, ntr=' num2str(numel(rasters(idxPdc).tr_idx))])
% fill([1:250 fliplr(1:250)],...
%     [FRsmooth_STANDARD+SDsmooth_STANDARD fliplr(FRsmooth_STANDARD-SDsmooth_STANDARD)],...
%     'k','FaceAlpha',0.3,'EdgeColor','none')
% fill([1:250 fliplr(1:250)],...
%     [FRsmooth_2+SDsmooth_2 fliplr(FRsmooth_2-SDsmooth_2)],...
%     'r','FaceAlpha',0.45,'EdgeColor','none')
% fill([1:250 fliplr(1:250)],...
%     [FRsmooth_6+SDsmooth_6 fliplr(FRsmooth_6-SDsmooth_6)],...
%     'b','FaceAlpha',0.45,'EdgeColor','none')
% plot(FRsmooth_STANDARD,'k','LineWidth',plotOptions.lineWidth)
% plot(FRsmooth_2,'r','LineWidth',plotOptions.lineWidth)
% plot(FRsmooth_6,'b','LineWidth',plotOptions.lineWidth)
% 
% if 1==(nRsp*2-1)
%     set(gca,'xtick',0:50:250); xlabel('time (ms)'); ylabel('FR (sp/s)')
% else
%     set(gca,'xtick',[],'ytick',[])
% end
% 
% 
% %% Now plot jitters on periodic
% 
% 
% for ij = 1:numel(Jitters)
%     
%     % skip periodic
%     if ij==idxPdc
%         continue
%     end
%     
%     
%     %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%     clear x_j y_j spikes FRsmooth
%     [x_j,y_j,prev250FR,prev100FR] = standardPd_getspikes(rasters(ij));
%     %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%     
%     spikes = zeros(max(y_j),250);
%     for iy = 1:max(y_j)
%         spikes(iy,x_j(y_j==iy)) = 1;
%     end
%     
%     % Calculate smoothed PSTH
%     FRsmooth_Jit = smoothPSTH_v2(spikes);
%     SDsmooth_Jit = smooth_STD(spikes,binsize);
%     
%     % Plot rasters
% %     figure(hf(1))
% %     subplot(nRsp,2,ij); hold on
% %     title([char(rasters(ij).jitter) ', ntr=' num2str(numel(rasters(ij).tr_idx))])
% %     for it = 1:(max(y)+max(y_j))
% %         plot([0 250],[it it],'k', 'LineWidth', 0.5*plotOptions.rasterLineWidth)
% %     end
% %     plot(x,y,'k+','MarkerSize',plotOptions.markerSize, 'LineWidth', plotOptions.rasterLineWidth)
% %     plot(x_j,y_j+max(y),'+',...
% %         'MarkerSize',plotOptions.markerSize,...
% %         'LineWidth', plotOptions.rasterLineWidth,...
% %         'Color',ALLcolors(7,:)); %strcmp(strtok(rasters(ij).jitter,'_'),strtok(plotOptions.colSelect,'_')), : ))
% %     ylim([0 max(y_j)+max(y)+1])
% %     
% %     
% %     % Plot PSTH - mean
% %     figure(hf(2))
% %     subplot(nRsp,2,ij); hold on
% %     plot(FRsmooth_STANDARD,'k','LineWidth',plotOptions.lineWidth)
% %     plot(FRsmooth_Jit,'LineWidth',plotOptions.lineWidth,...
% %         'Color',ALLcolors(7,:)); %strcmp(strtok(rasters(ij).jitter,'_'),strtok(plotOptions.colSelect,'_')), : ))
% %     title([char(rasters(ij).jitter) ', ntr=' num2str(numel(rasters(ij).tr_idx))])
% %     
% 
%     % Plot PSTH - with STD
%     figure(hf)
%     subplot(nRsp,2,ij); hold on
%     
%     fill([1:250 fliplr(1:250)],...
%         [FRsmooth_STANDARD+SDsmooth_STANDARD fliplr(FRsmooth_STANDARD-SDsmooth_STANDARD)],...
%         'k','FaceAlpha',0.3,'EdgeColor','none')
%     fill([1:250 fliplr(1:250)],...
%         [FRsmooth_Jit+SDsmooth_Jit fliplr(FRsmooth_Jit-SDsmooth_Jit)],...
%         ALLcolors(7,:),'FaceAlpha',0.45,'EdgeColor','none')
%     
%     plot(FRsmooth_STANDARD,'k','LineWidth',plotOptions.lineWidth)
%     plot(FRsmooth_Jit,'LineWidth',plotOptions.lineWidth,...
%         'Color',ALLcolors(7,:)); %strcmp(strtok(rasters(ij).jitter,'_'),strtok(plotOptions.colSelect,'_')), : ))
%     
%     title([char(rasters(ij).jitter) ', ntr=' num2str(numel(rasters(ij).tr_idx))])
%     
%     if ij==(nRsp*2-1)
%         set(gca,'xtick',0:50:250); xlabel('time (ms)'); ylabel('FR (sp/s)')
%     else
%         set(gca,'xtick',[],'ytick',[])
%     end
%     
%     
% end  %ij (jitters)
% 
% 
% % Set axis limits to be same for all subplots
% hAllAxes = findobj(hf,'type','axes');
% ylims = [min(cellfun(@min,get(hAllAxes,'YLim'))) max(cellfun(@max,get(hAllAxes,'YLim')))];
% set(hAllAxes,'YLim',ylims)
% set(hAllAxes,'XLim',[0 250])
% 
% 
% titlestr = sprintf('%s %s\n%i-%i Hz | %i dBSPL | BehavState: %s | %s',session,cluname,stimpars(1),stimpars(2),stimpars(3),behav,indVar);
% figure(hf); suptitle(titlestr)
% 
% savename = sprintf('%s_%s_%i_%i_%i_%i_%s_%s',session,cluname,stimpars(1),stimpars(2),stimpars(3),stimpars(4),behav,indVar);
% 
% savedir = fullfile(fn.processed,'PSTH_screening');
% if ~exist(savedir,'dir')
%     mkdir(savedir)
% end
% 
% set(gcf,'PaperPositionMode','auto','PaperOrientation','portrait')
% print(hf,fullfile(savedir,savename),'-depsc')
% 



end






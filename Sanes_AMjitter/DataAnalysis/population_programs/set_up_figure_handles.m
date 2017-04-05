function set_up_figure_handles(parname)

global Figs

% Only execute if the handles have not yet been created
if ~isfield(Figs,parname)
    
    % * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    % Neurometrics
    % * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    Figs.(parname).nmF = ['hF' num2str(numel(fieldnames(Figs))+1)];
    eval(sprintf('%s = figure;',Figs.(parname).nmF))
    
    scrsz = get(0,'ScreenSize');
    eval(sprintf('set(%s,''Position'',[1 scrsz(4) scrsz(3) scrsz(4)],''Nextplot'',''add'');',Figs.(parname).nmF))
    
    nsubplots = 1;
    for isp = 1:nsubplots
        Figs.(parname).nmS = ['hS' num2str(numel(fieldnames(Figs)))];
        eval(sprintf('%s=subplot(1,%i,%i);',Figs.(parname).nmS,nsubplots,isp))
    end
    
    
    % * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    % Histograms
    % * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

%     fighandle = ['hF' num2str(numel(fieldnames(Figs))+1) 'hist'];
%     Figs.(parname).histF = fighandle;
%     
%     eval(sprintf('%s = figure;',fighandle))
%     scrsz = get(0,'ScreenSize');
%     eval(sprintf('set(%s,''Position'',[1 scrsz(4) scrsz(3) scrsz(4)],''Nextplot'',''add'');',fighandle))
%     
%     switch indVar
%         case 'depth'
%             nsubplots = 4;
%         case 'jitter'
%             nsubplots = 10;
%     end
%     for isp = 1:nsubplots
%         SPhandle = ['hS' num2str(numel(fieldnames(Figs))) 'hist'];
%         Figs.(parname).histS = SPhandle;
%         eval(sprintf('%s=subplot(%i,1,%i);',SPhandle,nsubplots,isp))
%     end
%     
%     SPhandle = ['hS' num2str(numel(fieldnames(Figs))) 'nm'];
%     Figs.(parname).nmS = SPhandle;
%     
    
end


end





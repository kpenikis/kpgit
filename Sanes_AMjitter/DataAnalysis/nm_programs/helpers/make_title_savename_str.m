function [title_str,savename] = make_title_savename_str(parstruct,bs,subject,session,channel,clu,METRIC,binsize,indVar,method)

% Get pars for title and savename
HP    = parstruct(1).HP;
LP    = parstruct(1).LP;
dB    = parstruct(1).dB;
rate  = parstruct(1).AMrate;
block = parstruct(1).block;

% Add stimulus info in title
switch method
    
    case 'classif'
        title_str = sprintf('ch %i clu%i:  neurometric %s - %s bin %i  |  %s  |  noise: %i - %i Hz  |  %idB  |  AM %i Hz  |  blk%i',...
            channel,clu,indVar,METRIC,binsize,bs,HP,LP,dB,rate,block);
        
        savename  = sprintf('%s_%s_nm%s_classif_%s-bin%i_ch%i_clu%i_AM%iHz_%idB_%i-%i_blk%i',...
            subject,session,indVar,METRIC,binsize,channel,clu,rate,dB,HP,LP,block);
        
    case 'formulaFR'
        title_str = sprintf('ch %i clu%i:  neurometric %s - FR formula  |  %s  |  noise: %i - %i Hz  |  %idB  |  AM %i Hz  |  blk%i',...
            channel,clu,indVar,bs,HP,LP,dB,rate,block);
        
        savename  = sprintf('%s_%s_nm%s_formulaFR_ch%i_clu%i_AM%iHz_%idB_%i-%i_blk%i',...
            subject,session,indVar,channel,clu,rate,dB,HP,LP,block);
        
end

end

% Vary AM rates away from integer relationship
% 
%   2018-09-11


A = 0.6;

old_rates = 2.^[1 2 3 4 5];

new_rates = 2 .^ ([1 2 3 4 5] + A.*rand(1,5)-A/2);
new_pds   = round(1000./new_rates);
new_rates = 1000./new_pds;

figure; 
plot(sort([old_rates old_rates old_rates]),repmat([0 1 0],1,5),'k')
ax=gca;
ax.XScale = 'log';
hold on
plot(sort([new_rates new_rates new_rates]),repmat([0 1 0],1,5),'b')


pn = '/Users/kpenikis/Documents/SanesLab/LocalData/AM_IR_aversive/Stimuli';
savedir = fullfile(pn,'nonRatio');
if ~isdir(savedir)
    mkdir(savedir);
end

stimfiles = dir(fullfile(pn,'*.mat'));

for ii = 1:numel(stimfiles)
    
    load(fullfile(pn,stimfiles(ii).name))
    buffer_old = buffer;
    
    [~, i_mem] = ismember(buffer,old_rates);
    buffer = new_rates(i_mem);
    
    duration = sum(new_pds(i_mem));
    
    savename = strtok(stimfiles(ii).name,'.');
    save(fullfile(savedir,[strtok(stimfiles(ii).name,'.') '_' num2str(duration)]),'buffer','-v7.3')
    
    clear buffer buffer_old
end



% Also create new Warn stimulus with fast rate, to make duration easier to
% curate in circuit

NewWarnRate = 128;

Durations = [1000 1500 2000];
for duration = Durations
    
    buffer = repmat(NewWarnRate,1,1+duration/(1000/NewWarnRate));
    
    save(fullfile(pn,['rateVec_fastWarn' num2str(duration)]),'buffer','-v7.3')
    
end




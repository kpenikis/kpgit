

subject = 'WWWf_244303';
session = 'OA';


% Load data files
fn = set_paths_directories(subject,session);
fprintf('Session %s:  loading data...\n',session)
filename = sprintf( '%s_sess-%s_Spikes'   ,subject,session); load(fullfile(fn.processed,subject,filename));
filename = sprintf( '%s_sess-%s_Info'     ,subject,session); load(fullfile(fn.processed,subject,filename));
filename = sprintf( '%s_sess-%s_SoundData',subject,session); load(fullfile(fn.processed,subject,filename));
load(fullfile(fn.stim,'IRsequences.mat'))


%%

% Get block onset times

bkStarts_samp = find(diff(SoundData(8,:))~=0)+1;
bkStops_samp  = bkStarts_samp-1;
bkStarts_samp = bkStarts_samp(1:end-1);
bkStops_samp  = bkStops_samp(2:end);

b1=26;
b2=33;


% Plot it
AMrates = [2 4 8 16 32 64];
colors = [ 84  24  69;...
           120  10  41;...
           181   0  52;...
           255  87  51;...
           255 153   0;...
           255 205  60 ]./255;
colors = [colors; 0.7.*bone(4)];

set(0,'DefaultAxesFontSize',24)

scrsz = get(0,'ScreenSize');
figsize1 = [1 scrsz(4)/2 scrsz(3) scrsz(4)/2];

figure; hold on
set(gcf,'Position',figsize1)

hs1=subplot(6,1,1:3); hold on
set(gca,'yscale','log','xlim',[bkStarts_samp(b1) bkStops_samp(b2)])
box off
plot(bkStarts_samp(b1):bkStops_samp(b2),SoundData(1,bkStarts_samp(b1):bkStops_samp(b2)),...
    'k','LineWidth',2)
set(gca,'xtick',[],'ytick',[2 4 8 16 32 64])
ylabel('AM rate (Hz)')

hs2=subplot(6,1,4:6); hold on
set(gca,'yscale','linear','xlim',[bkStarts_samp(b1) bkStops_samp(b2)])
box off
set(gca,'xtick',[],'ytick',[])

for ib = b1:b2
        
    PdSwitches = find(diff(SoundData(1,bkStarts_samp(ib):bkStops_samp(ib)))~=0)+1;
    
    % periodic block
    if isempty(PdSwitches)
        
        subplot(hs1)
        plot(bkStarts_samp(ib):bkStops_samp(ib),SoundData(1,bkStarts_samp(ib):bkStops_samp(ib)),...
            'Color',colors(AMrates==SoundData(1,bkStarts_samp(ib)),:),'LineWidth',8)
        
        subplot(hs2)
        plot(bkStarts_samp(ib):bkStops_samp(ib),SoundData(2,bkStarts_samp(ib):bkStops_samp(ib)),...
            'Color',colors(AMrates==SoundData(1,bkStarts_samp(ib)),:))
    
    % IR block
    else  
        PdEdges = [bkStarts_samp(ib) bkStarts_samp(ib)+PdSwitches-1 bkStops_samp(ib)];
        for ir = 2:numel(PdEdges)
            
            subplot(hs1)
            plot(PdEdges(ir-1):(PdEdges(ir)-1),SoundData(1,PdEdges(ir-1):(PdEdges(ir)-1)),...
                'Color',colors(AMrates==SoundData(1,PdEdges(ir-1)),:),'LineWidth',8)
            
            subplot(hs2)
            plot(PdEdges(ir-1):(PdEdges(ir)-1),SoundData(2,PdEdges(ir-1):(PdEdges(ir)-1)),...
                'Color',colors(AMrates==SoundData(1,PdEdges(ir-1)),:))
            
        end        
    end
    
end












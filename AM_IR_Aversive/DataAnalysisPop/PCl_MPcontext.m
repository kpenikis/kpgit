function Data = PCl_MPcontext
%
%  PCl_MPcontext
%   Population Classifier of MP context
%
%
%  KP, 2019-09
%

% NEXT: 
%   separate match MP program from run classifier program
%   exclude ITI trials
%   re-run classifier for individual units now that ntrs is limited


global fn AMrates trMin Iterations

%!!!!!!!!!!!!!!!!!!!
Iterations  =  5000;
%!!!!!!!!!!!!!!!!!!!


%% Load Unit data files

fn = set_paths_directories;

q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q
%-------
spkshift = mean([UnitData([UnitData.IntTime_spk]>0).IntTime_spk]);
%-------


savedir = fullfile(fn.processed,'MPHclassifier');
q = load(fullfile(savedir,'ClassData'));
Data = q.Data;
clear q

savedir = fullfile(fn.processed,'MPHclassifier','Pop');
if ~exist(savedir,'dir')
    mkdir(savedir)
end


%% Which units to collect data from?

theseUnits = 1:size(Data,1);


%% Go through each MP comparison for each AM rate

Pop_dPrimes = nan(4,5);
mintrs      = nan(1,5);
padTrs      = 300;

for irate = 1:5
    
    fprintf('Getting data for %i Hz...\n',AMrates(irate))
    
    PdcRaster_7_1 = nan(padTrs,0);
    PdcRaster_7_2 = nan(padTrs,0);
    PdcRaster_8_1 = nan(padTrs,0);
    PdcRaster_8_2 = nan(padTrs,0);
    
    IrrRaster_7_1 = nan(padTrs,0);
    IrrRaster_7_2 = nan(padTrs,0);
    IrrRaster_8_1 = nan(padTrs,0);
    IrrRaster_8_2 = nan(padTrs,0);
    
    
    % Step through Units and collect all spikes
    
    h = waitbar(0,sprintf('Concatenating raster data from all units for %i Hz MPs',AMrates(irate)));
    
    for iUn = theseUnits
        
        % Determine which Irr MPs were presented to this unit
        % now it's a field
        if ~isfield(Data(iUn,irate).data,'IRid')
            fprintf('skipping iUn %i\n',iUn)
            continue
        end
        
        waitbar(iUn/numel(theseUnits))
        
        theseIrrMPs = [[Data(iUn,irate).data(:,1).IRid]; [Data(iUn,irate).data(:,1).iseq]]';
        
        for ic = 1:size(theseIrrMPs,1)
            
            % Skip if not enough trials (in either Pdc or Irr MP)
            if size(Data(iUn,irate).data(ic,1).raster,1)<trMin*2 || size(Data(iUn,irate).data(ic,2).raster,1)<trMin*2
                continue
            end
            
            % Concatenate Pdc raster data from this cell
            raster = nan(padTrs,ceil(1000/AMrates(irate)));
            raster(1:size(Data(iUn,irate).data(ic,1).raster,1),:) = Data(iUn,irate).data(ic,1).raster;
            
            eval(sprintf( 'PdcRaster_%i_%i = [PdcRaster_%i_%i raster];',theseIrrMPs(ic,1),theseIrrMPs(ic,2),theseIrrMPs(ic,1),theseIrrMPs(ic,2))) 
            
            
            % Concatenate Irr raster data from this cell
            raster = nan(padTrs,ceil(1000/AMrates(irate)));
            raster(1:size(Data(iUn,irate).data(ic,2).raster,1),:) = Data(iUn,irate).data(ic,2).raster;
            
            eval(sprintf( 'IrrRaster_%i_%i = [IrrRaster_%i_%i raster];',theseIrrMPs(ic,1),theseIrrMPs(ic,2),theseIrrMPs(ic,1),theseIrrMPs(ic,2))) 
            
        end %ic
        
    end %iUn
    
    close(h);
    
    % Trim NaNs
    PdcRaster_7_1 = PdcRaster_7_1(1:sum(~isnan(mean(PdcRaster_7_1,2))),:);
    PdcRaster_7_2 = PdcRaster_7_2(1:sum(~isnan(mean(PdcRaster_7_2,2))),:);
    PdcRaster_8_1 = PdcRaster_8_1(1:sum(~isnan(mean(PdcRaster_8_1,2))),:);
    PdcRaster_8_2 = PdcRaster_8_2(1:sum(~isnan(mean(PdcRaster_8_2,2))),:);
    IrrRaster_7_1 = IrrRaster_7_1(1:sum(~isnan(mean(IrrRaster_7_1,2))),:);
    IrrRaster_7_2 = IrrRaster_7_2(1:sum(~isnan(mean(IrrRaster_7_2,2))),:);
    IrrRaster_8_1 = IrrRaster_8_1(1:sum(~isnan(mean(IrrRaster_8_1,2))),:);
    IrrRaster_8_2 = IrrRaster_8_2(1:sum(~isnan(mean(IrrRaster_8_2,2))),:);
    
    % Get min n trials
    mintrs(irate) = min([size(PdcRaster_7_1,1) size(PdcRaster_8_1,1) size(IrrRaster_7_1,1) size(IrrRaster_8_1,1)]);
    
    % Put population data into struct
    clear PopMPH
    PopMPH = struct('Context',{'Pdc' 'IR'; 'Pdc' 'IR'; 'Pdc' 'IR'; 'Pdc' 'IR'},...
        'IRid',{7 7; 7 7; 8 8; 8 8},...
        'iseq',{1 1; 2 2; 1 1; 2 2},...
        'nTrs',{size(PdcRaster_7_1,1)  size(IrrRaster_7_1,1); ...
                size(PdcRaster_7_2,1)  size(IrrRaster_7_2,1); ...
                size(PdcRaster_8_1,1)  size(IrrRaster_8_1,1); ...
                size(PdcRaster_8_2,1)  size(IrrRaster_8_2,1) },...
        'raster',{PdcRaster_7_1 IrrRaster_7_1;...
                  PdcRaster_7_2 IrrRaster_7_2;...
                  PdcRaster_8_1 IrrRaster_8_1;...
                  PdcRaster_8_2 IrrRaster_8_2});
        %First draft, clip rasters to min N trials (later randomized which selected on each iteration)
    
    % Run classifier
    Pop_Res_L1o  = get_classifier_data( PopMPH, -1 );
    
    Pop_dPrimes(:,irate) = Pop_Res_L1o.dprime(:,2);
    
end %irate

save(fullfile(savedir,'Pop_dPrimes'),'Pop_dPrimes','-v7.3')



%% Compare population decoding to individual units

% Figure settings
set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',16)
scrsz = get(0,'ScreenSize');   %[left bottom width height]
fullscreen = [1 scrsz(4) scrsz(3) scrsz(4)];

% Each rate separately
hf=figure;
hold on
set(hf,'Position',fullscreen)
all_dPrimes = [];

for irate = 1:size(Data,2)
    each_dPrime = [];
    for iUn = 1:size(Data,1)
        if isempty(Data(iUn,irate).data) || ~isfield(Data(iUn,irate).data,'IRid')
            continue
        end
        for ip = 1:size(Data(iUn,irate).data,1)
            if size(Data(iUn,irate).data(ip,1).raster,1) < mintrs(irate) || size(Data(iUn,irate).data(ip,2).raster,1) < mintrs(irate)
                continue
            end
            each_dPrime = [each_dPrime; Data(iUn,irate).Res_L1o.dprime(ip,2)];
            all_dPrimes = [all_dPrimes; Data(iUn,irate).Res_L1o.dprime(ip,2)];
        end %ip
    end %iUn
    subplot(1,5,irate); hold on
    plot([0 length(each_dPrime)+1],[0 0],':k')
    plot(sort(each_dPrime),'.k')
    plot((1:numel(Pop_dPrimes(:,irate)))+length(each_dPrime)/2,Pop_dPrimes(:,irate),'ob','MarkerSize',10)
    plot([0 length(each_dPrime)+1],[mean(Pop_dPrimes(:,irate)) mean(Pop_dPrimes(:,irate))],'-b','LineWidth',2)
    axis square
    xlim([0 length(each_dPrime)+1])
    ylim([-3 4])
    title([num2str(AMrates(irate)) ' Hz'])
    if irate==1
        xlabel('Pdc-Irr MP comarison #')
        ylabel('dprime')
    end
end %irate

print_eps_kp(hf,fullfile(fn.figs,'MPHclass','MPcontext_all-vs-pop'))

% All rates together
hfa=figure;
plot([0 length(all_dPrimes)+1],[0 0],':k')
hold on
plot(sort(all_dPrimes),'.k')
plot((1:numel(Pop_dPrimes))+length(all_dPrimes)/2,Pop_dPrimes(:),'ob','MarkerSize',10)
plot([0 length(all_dPrimes)+1],[mean(Pop_dPrimes(:)) mean(Pop_dPrimes(:))],'-b','LineWidth',2)
axis square
xlim([0 length(all_dPrimes)+1])
ylim([-3 4])
title('All rates and MPs')
xlabel('MP comarison #')
ylabel('dprime')

print_eps_kp(hfa,fullfile(fn.figs,'MPHclass','MPcontext_all-vs-pop_allrates'))

keyboard


end %function





function Wave = ap_stimplotting(subject,raster)
%  
global fn

% Load first block rate vector stimulus info
block = raster(1).block;
block = block(1);
stfs = dir(fullfile(fn.raw,subject,sprintf('Block-%i_Stim',block),'*.mat'));
stfns = {stfs.name};

% Get vectors of rates for stimuli of this block
for iif = 1:numel(stfs)
    rateVecs(iif) = load(fullfile(stfs(iif).folder,stfns{iif}));
end

Wave=struct();
for ks = 1:numel(raster)
    
    data = raster(ks);
    data.block = data.block(1);
    
    if block ~= data.block
        
        % Load new rate vector stimulus info
        block = data.block;
        stfs = dir(fullfile(stimdir,subject,sprintf('Block-%i_Stim',block),'*.mat'));
        stfns = {stfs.name};
        
        % Get vectors of rates for stimuli of this block
        for iif = 1:numel(stfs)
            rateVecs(iif) = load(fullfile(stfs(iif).folder,stfns{iif}));
        end
        
    end
    
    % Get rate vector for current stimulus
    rV = rateVecs(strcmp(stfns,data.stimfn)).buffer;
    rV = rV(2:end);
    
    fs = 1e3;
    t  = cumsum([0 1./rV]);
    s = round(t.*fs);
    y  = [];
    for ip = 1:numel(rV)
        tp=[]; tp = 0:(s(ip+1)-s(ip)-1);
        y = [y 1- data.AMdepth*cos(2*pi*rV(ip)/fs*tp)];
    end
%     hold on; plot(y)

    % Remove first 1/4 of the first period
    yclip = y(round(0.25/rV(1)*fs):end);
    
    % Add sound during unmodulated portion
    y_um = ones(1,round(data.AMonset/1000*fs)); %0.8165.*
    stim = [y_um yclip];
    
    % Plot individual stim, for debugging
%     hF(ks) = figure;
%     scrsz = get(0,'ScreenSize');
%     set(hF(ks),'Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)/2],...
%         'Nextplot','add');
%     fill([1:length(stim) length(stim):-1:1] /fs*1000,[stim fliplr(-stim)],[0.3 0.3 0.3])
%     hold on
%     xtx = 500:500:2000;
%     for itx=xtx
%         plot([itx itx],[-2 2],'k--')
%     end
%     xlim([0 length(stim)]/fs*1000)
%     xlabel('time (ms)')
%     title(data.stimfn)
    

    % Print comparison between Wav and Stim durations, as a check
%     fprintf('\nraster %i \nWav dur %i, Stim dur %i \nWav AM %i, Stim AM %i',...
%         ks,(length(stim)/fs*1000), data.stimDur, length(y_um), data.AMonset)
    

    % Save data to output struct
    Wave(ks).x_ms = [0 length(stim)]/fs*1000;
    Wave(ks).y    = stim;
    Wave(ks).fs   = fs;
    Wave(ks).stfn = data.stimfn;
    
end

end





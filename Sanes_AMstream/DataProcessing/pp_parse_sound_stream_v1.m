function SoundData2 = pp_parse_sound_stream_v1(SoundData,Info)
%
% need several outputs: 
%  - SoundData sampled at 1 ms
%      with block type label (and silence and unmod)


global fn

stimfolder = strsplit(Info.stimfiles, '\');
stimfolder = stimfolder{end};
stimfolder = fullfile(fn.stim,stimfolder);

stimfiles = dir(fullfile(stimfolder,'*.mat'));

if numel(stimfiles)==2
    protocol = 'Linearity';
elseif numel(stimfiles)==3
    protocol = 'SpectralSwitch';
elseif numel(stimfiles)==311
    protocol = 'Trials';
end

switch protocol
    
    case 'Linearity'
        q=load(fullfile(stimfolder,'fullStream_AMrateVec_26trs_20160323.mat'));
        pd_AMrates = q.buffer;
        clear q
        
        q=load(fullfile(stimfolder,'fullStream_BlockVec_26trs_20160323.mat'));
        pd_Blocks  = q.buffer;
        blockKey   = q.stim_array;
        clear q
        
        
        %% Make 1 kHz sampled versions of these period vectors
        
        % First set some parameters
        periodic_rates = [2 4 8 16 32 64];
        IR_dur = 2*sum(1./periodic_rates);
        periodic_nCycles = ceil(periodic_rates*IR_dur);
        pdc_durs = periodic_nCycles .* (1./periodic_rates);
        
        align_fs = 1000;%Info.fs_sound;
        
        % Create empty vectors for output
        ms_AMrates=[];
        ms_Blocks =[];
        
        % Get block transitions and go through each block
        newBlock = [1 1+find(diff(pd_Blocks))];
        for ib = 1:(numel(newBlock)-1)
            
            theseRates = pd_AMrates(newBlock(ib):(newBlock(ib+1)-1));
            thisBlock  = pd_Blocks(newBlock(ib));
            
            if thisBlock>=1 && thisBlock<=6
                
                ms_AMrates = [ms_AMrates repmat(theseRates(1),1,round(IR_dur*align_fs))];
                ms_Blocks  = [ms_Blocks  repmat(thisBlock,1,round(IR_dur*align_fs))];
                
            else
                for ir = theseRates
                    add_this=[];
                    add_this = repmat(ir,1,round((1000/ir)+(rand(1,1)/1000-0.0005)));
                    ms_AMrates = [ms_AMrates add_this];
                    ms_Blocks  = [ms_Blocks  repmat(thisBlock,1,length(add_this))];
                end
            end
                        
        end
        
        aaa=234;
        
        
        %% Make 1 kHz sampled versions of SoundData
        
        % Down sample SoundData 
        SoundData2 = SoundData(:,1+find(diff(floor([1:size(SoundData,2)]/Info.fs_sound*1000))));
        
        % Make sure matches sampling rate
        if abs( (size(SoundData,2)/size(SoundData2,2)) - round(Info.fs_sound)/1000 ) > 10^-3
            keyboard
        end
        
        id=0; comb=100;
        for ii = 1000:comb:(size(SoundData2,2)-size(ms_AMrates,2))
            id=id+1;
%             vecDist(id) = sum(SoundData2(1,ii:(ii+size(ms_AMrates,2)-1)) - ms_AMrates);
            vecDist(id) = sum(SoundData2(1,ii:(ii+2999)) - ms_AMrates(1,1:3000));
        end
        
        [~,im] = min(abs(vecDist)); im_ms = 1000+((im-1)*comb);
        
        comb=500;
        id=0; 
        vecDist2=nan(1,(im_ms+comb)-(im_ms-comb));
        for ii = (im_ms-comb):(im_ms+comb)
            id=id+1;
            vecDist2(id) = sum(SoundData2(1,ii:(ii+4999)) - ms_AMrates(1,1:5000));
        end
        
        [~,im2] = min(abs(vecDist2)); im_ms2 = im2 + im_ms-comb-1;
        
        
        if ~any(find(vecDist2==0))
            keyboard
        end
        
        idxALIGN = find(vecDist2==0) + (im_ms-comb) - 1;
        
        sum(SoundData2(1,idxALIGN:(idxALIGN+size(ms_AMrates,2)-1)) - single(ms_AMrates))
        
        % SUM IS STILL LARGE
        % is there a problem with subtracting double from single precision
        % number? no..
        % most likely that the # of ms for periods (see 16 especially)
        % aren't the same
        
        
        % Align them....
        
        
    case 'SpectralSwitch'
    case 'Trials'
end






end







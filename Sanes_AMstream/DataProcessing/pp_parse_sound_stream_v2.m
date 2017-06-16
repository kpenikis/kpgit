function [SoundData,Info] = pp_parse_sound_stream_v2(SoundData,Info)
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
        
        
        %% Make X kHz sampled versions of these period vectors
        
        % First set some parameters
        periodic_rates = [2 4 8 16 32 64];
        IR_dur = 2*sum(1./periodic_rates);
        periodic_nCycles = ceil(periodic_rates*IR_dur);
        pdc_durs = periodic_nCycles .* (1./periodic_rates);
        
        align_fs = Info.fs_sound;
        
        % Create empty vectors for output
        ms_AMrates=[];
        ms_Blocks =[];
        
        % Get block transitions and go through each block
        newBlock = [1 1+find(diff(pd_Blocks))];
        for ib = 1:(numel(newBlock)-1)
            
            theseRates = pd_AMrates(newBlock(ib):(newBlock(ib+1)-1));
            thisBlock  = pd_Blocks(newBlock(ib));
            
            if thisBlock>=1 && thisBlock<=6
                
                ms_AMrates = [ms_AMrates repmat(theseRates(1),1,ceil(IR_dur*align_fs))];
                ms_Blocks  = [ms_Blocks  repmat(thisBlock,1,ceil(IR_dur*align_fs))];
                
            else
                for ir = theseRates
                    add_this=[];
                    add_this = repmat(ir,1,ceil(align_fs/ir));
                    ms_AMrates = [ms_AMrates add_this];
                    ms_Blocks  = [ms_Blocks  repmat(thisBlock,1,length(add_this))];
                end
            end
                        
        end
        
        
        
        %%  Get to the next step.
        
        % X from saved AM rates, label block (refer to pdBlocks)
        %   from rms of sound, label unmodulated portion
        %   from rms of sound, label silence
        %   *maybe* change fs to 1 kHz at end
        
                
        AMrateTrans_input = [find(diff(ms_AMrates)) size(ms_AMrates,2)] +1;
        
        AMrateTrans_output = find(diff(SoundData(1,:)))+1;
        Blocks_SoundData = zeros(1,size(SoundData,2));
        
        try
        for iBk = 3:numel(AMrateTrans_output)
            
            if (iBk-1) > size(AMrateTrans_input,2)
                break
            end
            
            rates_out = unique(SoundData(1,AMrateTrans_output(iBk-1):(AMrateTrans_output(iBk)-1)),'stable');
            samples_out = numel(SoundData(1,AMrateTrans_output(iBk-1):(AMrateTrans_output(iBk)-1)));
            
            rates_in  = unique(ms_AMrates(1,AMrateTrans_input(iBk-2):(AMrateTrans_input(iBk-1)-1)),'stable');
            samples_in  = numel(ms_AMrates(1,AMrateTrans_input(iBk-2):(AMrateTrans_input(iBk-1)-1)));
            
            if rates_out ~= single(rates_in)
                sprintf('input rate: %i, recorded rate: %i',rates_in,rates_out)
                keyboard
            end
            
            block_in  = unique(ms_Blocks(1,AMrateTrans_input(iBk-2):(AMrateTrans_input(iBk-1)-1)),'stable');
            
            if numel(block_in)==1
                Blocks_SoundData(1,AMrateTrans_output(iBk-1):(AMrateTrans_output(iBk)-1)) = repmat(block_in,1,samples_out);
                
            elseif numel(block_in)==2
                
                sampsBk1 = round(samples_out*sum(ms_Blocks(1,AMrateTrans_input(iBk-2):(AMrateTrans_input(iBk-1)-1))==block_in(1))/samples_in);
                
                Blocks_SoundData(1,[1:sampsBk1]+AMrateTrans_output(iBk-1)-1) = repmat(block_in(1),1,sampsBk1);
                Blocks_SoundData(1,[(1+sampsBk1):samples_out]+AMrateTrans_output(iBk-1)-1) = repmat(block_in(2),1,samples_out-sampsBk1);
                
            else
                warning('too many blocks to handle right now')
                keyboard
            end
            
        end
        
        catch
            keyboard
        end
        
        
        
        %%   from rms of sound, label unmodulated portion
        sound_out = SoundData(2,1:AMrateTrans_output(2));
        
        
        rmsbinsize = round(0.2*Info.fs_sound);
        for ii = 1:AMrateTrans_output(2)
            sound_out(ii) = rms(SoundData(2,ii:(ii+rmsbinsize)));
        end
        
        coeffFilt   = ones(1,rmsbinsize)/rmsbinsize;
        
        sound_smooth = filter(coeffFilt, 1, sound_out );
        
        fDelay    = (length(coeffFilt)-1)/2;
        adjustedS = [1:length(sound_smooth)]-fDelay;
                
        [pks,ipks] = findpeaks(double(diff(sound_smooth)./max(diff(sound_smooth))),'SortStr','descend','MinPeakProminence',0.3);
        
        unmON = find(SoundData(2,1:AMrateTrans_output(2)) > 0);
        unmON = unmON(1);
        
        %
%         figure(100); clf; hold on
%         plot(adjustedS,sound_out./max(sound_out),'LineWidth',2)
%         plot(adjustedS(2:end),diff(sound_smooth)./max(diff(sound_smooth)),'r','LineWidth',2)
%         plot(adjustedS(ipks(2)),pks(2),'*k','MarkerSize',10,'LineWidth',3)
        %
        hf=figure; hold on 
        plot(SoundData(2,1:AMrateTrans_output(2)),'Color',[0.4 0.4 0.4])
        plot([1 ipks(2)+1 ipks(2)+1 AMrateTrans_output(2)-1 AMrateTrans_output(2)-1 AMrateTrans_output(2)],[0 0 1 1 0 0],'--','LineWidth',2,'Color',[255,140,0]./255)
        plot([1 unmON unmON ipks(2) ipks(2) AMrateTrans_output(2)],[0 0 1 1 0 0],'--g','LineWidth',2)
        plot([1 AMrateTrans_output(1) AMrateTrans_output(1) unmON-1 unmON-1 AMrateTrans_output(2)],[0 0 1 1 0 0],'--k','LineWidth',2)
        %
        
        
        result = input('\ndoes the plot look ok?');
        if ~isempty(result)
            keyboard
        else 
            close(hf)
        end        
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % Add silent period before sound onset to SD block label
        Blocks_SoundData(1, AMrateTrans_output(1) : (unmON-1) ) = 12;
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~           
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % Add unmodulated sound (front end) to SD block label
        Blocks_SoundData(1, unmON : ipks(2) ) = 11;
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % Add first block code, now that timestamp is set
        Blocks_SoundData(1, (ipks(2)+1) : (AMrateTrans_output(2)-1) ) = find(strcmp(num2str(SoundData(1,AMrateTrans_output(1))),blockKey));
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        
        
        %% Unmod and silence at back end
        
                
        post_AM =  SoundData(2, AMrateTrans_output(end-1):end) ~= 0;
        silence_end = SoundData(2, AMrateTrans_output(end-1):end) == 0;
        
        %
        hf=figure; hold on
        plot( SoundData(2, AMrateTrans_output(end-1):end) ,'Color',[0.4 0.4 0.4])
        plot( SoundData(2, AMrateTrans_output(end-1):end) ~= 0 , '--g', 'LineWidth',2 ) %unmod
        plot( SoundData(2, AMrateTrans_output(end-1):end) == 0 , '--k', 'LineWidth',2 )%silence
        xlim([1 size(SoundData,2)-AMrateTrans_output(end-1)])
        %
        
        if sum(diff(post_AM)) ~= -1
            warning('unmod segment not necessarily clean')
            keyboard
        end
        if sum(diff(silence_end)) ~= 1
            warning('silent segment not necessarily clean')
            keyboard
        end
        
        result = input('\ndoes the plot look ok?\n');
        if ~isempty(result)
            keyboard
        else 
            close(hf)
        end
        
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % As long as the UNMODULATED portion is longer than 5 seconds, add
        % it to Blocks_SoundData
        if length(find(post_AM)) > (Info.fs_sound*5) 
            Blocks_SoundData(1, find(post_AM)+AMrateTrans_output(end-1)-1 ) = 11;
        end
        
        % As long as the SILENCE portion is longer than 5 seconds, add
        % it to Blocks_SoundData
        if length(find(post_AM)) > (Info.fs_sound*5) 
            Blocks_SoundData(1, find(silence_end)+AMrateTrans_output(end-1)-1 ) = 12;
        end
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        
        
        %% Add Blocks_SoundData to SoundData
        
        SoundData = [SoundData; Blocks_SoundData];
        
        
        %% Add unmodulated and silent codes to the block key, and save to Info struct
        
        if numel(blockKey) ~= 10
            keyboard
        else
            blockKey{11} = 'unmod';
            blockKey{12} = 'silence';
            Info.blockKey = blockKey;
        end
        
        Info.sound_rows{8} = 'Blocks';
        
        
        
        %% Make 1 kHz sampled versions of SoundData
        
        
    case 'SpectralSwitch'
        keyboard
        
    case 'Trials'
        keyboard
        
end




end







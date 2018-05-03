function [SoundData,Info] = get_block_labels(pd_AMrates,pd_Blocks,blockKey,SoundData,Info)
% Called by pp_parse_sound_stream_v2, which is called within pp_make_InfoPhys
% 

%% Make X kHz sampled versions of these period vectors
%   (match fs of SoundData)

% First set some parameters
periodic_rates = [2 4 8 16 32 64];
IR_dur = 2*sum(1./periodic_rates);
periodic_nCycles = ceil(periodic_rates*IR_dur);
pdc_durs = periodic_nCycles .* (1./periodic_rates);

align_fs = Info.fs_sound;

% Create empty vectors for output
samp_AMrates_in=[];
samp_Blocks_in =[];

% Get block transitions and go through each block
newBlock = [1 1+find(diff(pd_Blocks))];
for ib = 1:(numel(newBlock)-1)
    
    theseRates = pd_AMrates(newBlock(ib):(newBlock(ib+1)-1));
    thisBlock  = pd_Blocks(newBlock(ib));
    
    if thisBlock>=1 && thisBlock<=6
        
        if abs(length(theseRates)-periodic_nCycles(thisBlock))>1 && length(theseRates)/periodic_nCycles(thisBlock) ~= 2 
            %difference of n periods more than one but not double
            keyboard
        end
        
        samp_AMrates_in = [samp_AMrates_in repmat(theseRates(1),1,ceil(pdc_durs(thisBlock)*align_fs))];
        samp_Blocks_in  = [samp_Blocks_in  repmat(thisBlock,1,ceil(pdc_durs(thisBlock)*align_fs))];
        
        if length(theseRates)/periodic_nCycles(thisBlock) == 2
            samp_AMrates_in = [samp_AMrates_in repmat(theseRates(1),1,ceil(pdc_durs(thisBlock)*align_fs))];
            samp_Blocks_in  = [samp_Blocks_in  repmat(thisBlock,1,ceil(pdc_durs(thisBlock)*align_fs))];
        end
        
    else
        for ir = theseRates
            add_this=[];
            add_this = repmat(ir,1,ceil(align_fs/ir));
            samp_AMrates_in = [samp_AMrates_in add_this];
            samp_Blocks_in  = [samp_Blocks_in  repmat(thisBlock,1,length(add_this))];
        end
    end
    
end



%%
%   from saved AM rates, label block (refer to pdBlocks)
%   from rms of sound, label unmodulated portion
%   from rms of sound, label silence

AMrateTrans_input = [find(diff(samp_AMrates_in)) size(samp_AMrates_in,2)] +1;

AMrateTrans_output = find(diff(SoundData(1,:)))+1;
Blocks_SoundData = zeros(1,size(SoundData,2));

try
    for iBk = 3:numel(AMrateTrans_output) % go through each AM rate transition recorded in SoundData
        
        if (iBk-1) > size(AMrateTrans_input,2)
            break
        end
        
        rates_out = unique(SoundData(1,AMrateTrans_output(iBk-1):(AMrateTrans_output(iBk)-1)),'stable');
        samples_out = numel(SoundData(1,AMrateTrans_output(iBk-1):(AMrateTrans_output(iBk)-1)));
        
        rates_in  = unique(samp_AMrates_in(1,AMrateTrans_input(iBk-2):(AMrateTrans_input(iBk-1)-1)),'stable');
        samples_in  = numel(samp_AMrates_in(1,AMrateTrans_input(iBk-2):(AMrateTrans_input(iBk-1)-1)));
        
        if rates_out ~= single(rates_in)
            sprintf('input rate: %i, recorded rate: %i',rates_in,rates_out)
            keyboard
        end
        if abs(samples_out-samples_in)>50
            
            keyboard
            
%             round(samples_out/Info.fs_sound*1000)
% %             
%             figure; clf;hold on
%             plot(SoundData(1,AMrateTrans_output(iBk-1):(AMrateTrans_output(iBk)-1)),'r')
%             plot(samp_Blocks_in(1,AMrateTrans_input(iBk-2):(AMrateTrans_input(iBk-1)-1)),'b')
%             plot(samp_AMrates_in(1,AMrateTrans_input(iBk-2):(AMrateTrans_input(iBk-1)-1)),'k')
%             
%             figure;
%             plot(SoundData(1,AMrateTrans_output(iBk-5):(AMrateTrans_output(iBk+5)-1)),'g','LineWidth',2)
            
        end
        
        block_in  = unique(samp_Blocks_in(1,AMrateTrans_input(iBk-2):(AMrateTrans_input(iBk-1)-1)),'stable');
        
        if numel(block_in)==1
            Blocks_SoundData(1,AMrateTrans_output(iBk-1):(AMrateTrans_output(iBk)-1)) = repmat(block_in,1,samples_out);
            
        elseif numel(block_in)==2
            
            sampsBk1 = round(samples_out * sum(samp_Blocks_in(1,AMrateTrans_input(iBk-2):(AMrateTrans_input(iBk-1)-1))==block_in(1)) /samples_in);
            
            Blocks_SoundData(1,[1:sampsBk1] + AMrateTrans_output(iBk-1)-1) = repmat(block_in(1),1,sampsBk1);
            Blocks_SoundData(1,[(1+sampsBk1):samples_out] + AMrateTrans_output(iBk-1)-1) = repmat(block_in(2),1,samples_out-sampsBk1);
            
            if block_in(1)<block_in(2) && ( abs((sampsBk1/Info.fs_sound)-2)<0.1 || abs((sampsBk1/Info.fs_sound)-4)<0.1 )
                aaa='good';
            elseif block_in(2)<block_in(1) && abs(((samples_out-sampsBk1)/Info.fs_sound)-2)<0.1
                aaa='good';
            else
                keyboard
                figure;
                plot(Blocks_SoundData(1,[1:samples_out] + AMrateTrans_output(iBk-1)-1),'LineWidth',2)
                hold on
                plot(SoundData(2,[1:samples_out] + AMrateTrans_output(iBk-1)-1),'k','LineWidth',2)
            end
            
            
        else
            warning('too many blocks to handle right now')
            keyboard
        end
        
    end
    
catch
    keyboard
end



%%  Label unmodulated and silent periods in beginning of session
%     timepoints based on rms of sound

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


%         figure(100); clf; hold on
%         plot(adjustedS,sound_out./max(sound_out),'LineWidth',2)
%         plot(adjustedS(2:end),diff(sound_smooth)./max(diff(sound_smooth)),'r','LineWidth',2)
%         plot(adjustedS(ipks(2)),pks(2),'*k','MarkerSize',10,'LineWidth',3)
%         plot(adjustedS(unmON),pks(2),'*k','MarkerSize',10,'LineWidth',3)


hf=figure; hold on
plot(SoundData(2,1:AMrateTrans_output(2)),'Color',[0.4 0.4 0.4])
ip(1)=plot([1 AMrateTrans_output(1) AMrateTrans_output(1) unmON-1 unmON-1 AMrateTrans_output(2)],[0 0 1 1 0 0],'--k','LineWidth',2); %silence
ip(2)=plot([1 unmON unmON ipks(2) ipks(2) AMrateTrans_output(2)],[0 0 1 1 0 0],'--g','LineWidth',2); %unmodulated noise
ip(3)=plot([1 ipks(2)+1 ipks(2)+1 AMrateTrans_output(2)-1 AMrateTrans_output(2)-1 AMrateTrans_output(2)],[0 0 1 1 0 0],'--','LineWidth',2,'Color',[255,140,0]./255); %AM

legend(ip,{'silence' 'unmod.' 'AM'})

%


result = input('\ndoes the plot look ok?');
if ~isempty(result)
    keyboard
else
    close(hf)
    clear ip
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



%% Label unmodulated and silent periods at end of session


UM_end =  SoundData(2, AMrateTrans_output(end-1):end) ~= 0;
silence_end = SoundData(2, AMrateTrans_output(end-1):end) == 0;

%
hf=figure; hold on
plot( SoundData(2, AMrateTrans_output(end-1):end) ,'Color',[0.4 0.4 0.4])
ip(1)=plot( UM_end      , '--g', 'LineWidth',2 ); %unmod
ip(2)=plot( silence_end , '--k', 'LineWidth',2 ); %silence
%         plot( SoundData(2, AMrateTrans_output(end-2):end) ,'Color',[0.4 0.4 0.4])
xlim([1 size(SoundData,2)-AMrateTrans_output(end-1)])
legend(ip,{'unmod.' 'silence'})
%

%
% The current problem is that the block labeling could be buggy now that
% there are fewer unique blocks and more repeated transitions. Also not
% sure if unmodulated period at end is being labeled correctly now - it
% seems to start appx 1 block to early, perhaps because the last rate of an
% IR block is the same rate as the last periodic block.
% 

%
% hf=figure; hold on
% plot( SoundData(2, AMrateTrans_output(end-25):end) ,'Color',[0.4 0.4 0.4])
% plot( SoundData(1, AMrateTrans_output(end-25):end))
% plot( Blocks_SoundData(1, AMrateTrans_output(end-25):end))
%

if sum(diff(UM_end)) ~= -1
    warning('unmod segment not necessarily clean')
    keyboard
end
if sum(diff(silence_end)) ~= 1
    warning('silent segment not necessarily clean')
    keyboard
end

result = input('\ndoes the plot look ok?\ntype _0_ if ended early, or _4_ if unmod label starts too early\n');
if isempty(result)
    close(hf)
    clear ip
elseif result==4 || result==0
%     keyboard
else
    keyboard
end


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% As long as the UNMODULATED portion is longer than 1 second, add
% it to Blocks_SoundData
UM_end_samples = find(UM_end)+AMrateTrans_output(end-1)-1;
if       isempty(result)  &&  ( length(find(UM_end)) > (Info.fs_sound*1) ) 
    Blocks_SoundData(1, UM_end_samples ) = 11;

elseif   result==0  %session ended early
    Blocks_SoundData(1, UM_end_samples ) = 0;

elseif   result==4  %unmod starts later than marked
    
    delayed_UM_samples = UM_end_samples(2*Info.fs_sound:end);
    plot([delayed_UM_samples(1)-AMrateTrans_output(end-1)+1 delayed_UM_samples(1)-AMrateTrans_output(end-1)+1],[-5 5],'--c','LineWidth',2)
    check = input(' ...now does the unmod starttime look ok?');
    
    if isempty(check) && length(delayed_UM_samples) >= (Info.fs_sound*1)
        Blocks_SoundData(1, delayed_UM_samples ) = 11;
    elseif  isempty(check) && length(delayed_UM_samples) < (Info.fs_sound*1)
        Blocks_SoundData(1, delayed_UM_samples ) = 0;
    else
        keyboard
%         delayed_UM_samples = delayed_UM_samples(0.5*Info.fs_sound:end);
%         check=[];
    end
    close(hf)
    clear ip
end

% As long as the SILENCE portion is longer than 5 seconds, add
% it to Blocks_SoundData
if ( length(find(UM_end)) > (Info.fs_sound*1) ) && isempty(result)
    Blocks_SoundData(1, find(silence_end)+AMrateTrans_output(end-1)-1 ) = 12;
elseif result==4
    Blocks_SoundData(1, find(silence_end)+AMrateTrans_output(end-1)-1 ) = 0;
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

Info.sound_rows{end+1} = 'Blocks';


% figure; 
% plot(SoundData([1 2 8],:)','LineWidth',2)


end

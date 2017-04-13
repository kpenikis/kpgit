function drafting_irreg_AM_chunks
% This stimulus recipe allows for 3 important criteria to be matched:
%
%   1) a few specific AM rate vector sequences (chunks) are repeated for
%   multiple (15) trials, allowing for psth and trial-to-trial variability
%   analyses.
%
%   2) the specific chunks can arrive either early or late in an irregular
%   segment, allowing us to look for adaptation on a constant physical
%   signal.
%
%   3) each unique AM rate (/period duration) is presented in several
%   diverse positions within an irregular sequence. ie, there are several
%   unique AM rates  preceeding it. this allows for a rough assessment of
%   details behind any observed hysteresis effects, or conditions for
%   non-linearities, and allows for an assessment of neural activity
%   related to irregularity as a condition, with less bias by responses
%   specific to an individual sequence.
%
%   weaknesses:
%   - there may not be enough randomized chunks, so analyses of adaptation
%   may be confounded by responses specific to the rate sequences
%   presented.
%
%   each unique rate value is presented 120 times in the irregular
%   condition.
%

% rng('shuffle')

AMrates = 2.^[1:0.5:6];  %2-64 Hz
periodic_rates = [2 4 8 16 32];

irreg_chunk_dur = 2*sum(1./AMrates);

try
    load('example_irregular_sequence.mat')
    
catch
    SATISFIED = zeros(size(AMrates));
    while all(SATISFIED==0)
        
        randChunkA = randperm(length(AMrates));
        randChunkB = randperm(length(AMrates));
        randChunkC = randperm(length(AMrates));
        randChunkD = randperm(length(AMrates));
        rateVec_idx = [randChunkA randChunkB randChunkC randChunkD];
        
        while any(abs(diff(rateVec_idx))<2) %ensure there are no stepwise increments
            randChunkA = randperm(length(AMrates));
            randChunkB = randperm(length(AMrates));
            randChunkC = randperm(length(AMrates));
            randChunkD = randperm(length(AMrates));
            rateVec_idx = [randChunkA randChunkB randChunkC randChunkD];
        end
        
        for ir = 1:numel(AMrates) %ensure that each individual rate follows a variety of preceding rates
            idxs=[]; idxs = find(rateVec_idx==ir);
            if numel(unique([ ir rateVec_idx(idxs(idxs~=1)-1) ])) == (1+numel(idxs(idxs~=1)))
                SATISFIED(ir) = 1;
            end
        end
        if any(SATISFIED==0) %if not all do, start over
            SATISFIED = zeros(size(AMrates));
        end
    end
end

% Final AM rate "segments"
rateVec1    = AMrates([randChunkA randChunkB]);
rateVec1REV = AMrates([randChunkB randChunkA]);
rateVec2    = AMrates([randChunkC randChunkD]);
rateVec2REV = AMrates([randChunkD randChunkC]);

% save('example_irregular_sequence','AMrates','randChunkA','randChunkB','randChunkC','randChunkD','-v7.3')

% Plot the AM rate sequences
figure; hold on
xstart = 0;
plot(xstart+1:xstart+length(rateVec1),rateVec1,'+k')
xstart = xstart+length(rateVec1);
plot(xstart+1:xstart+length(rateVec1REV),rateVec1REV,'+b')
xstart = xstart+length(rateVec1REV);
plot(xstart+1:xstart+length(rateVec2),rateVec2,'+k')
xstart = xstart+length(rateVec2);
plot(xstart+1:xstart+length(rateVec2REV),rateVec2REV,'+b')
set(gca,'yscale','log')
set(gca,'ytick',AMrates)
xlabel('period number')
ylabel('AM rate (Hz)')



%% Estimating stimulus time, number of trials for a session

ml_consumed = 2.5;%ml
ml_flowrate = 0.25;%ml/m
total_time = ml_consumed / ml_flowrate * 60;

n_irreg_trials = 22;
irreg_repeats_time = irreg_chunk_dur * 4 * n_irreg_trials;

n_periodic_rates = numel(periodic_rates);
n_periodic_trials = 22;
periodic_time = n_periodic_rates * irreg_chunk_dur * n_periodic_trials;

standard_stim_time = irreg_repeats_time + periodic_time


time_left_for_rand_irregs = total_time - standard_stim_time;
n_rand_irregs = time_left_for_rand_irregs/irreg_chunk_dur;



%% Piece together long AM vector for continuous sound stream

periodic_rates;
rateVec1   ;
rateVec1REV;
rateVec2;
rateVec2REV;

REGblocks = 1:numel(periodic_rates);
IRblocks  = numel(periodic_rates)+(1:2);

stim_array = strsplit(num2str(periodic_rates));
stim_array{end+1} = 'rateVec_';
stim_array{end+1} = 'rateVec_REV';

%%

transitions = [];

for ir = REGblocks([2 4])  %allows for 11 trials of each type of transition
    for ii = IRblocks
        transitions = [transitions; ir ii];
    end
end





%%

this_random_pick = REGblocks(round(1 + (length(REGblocks)-1).*rand(1,1)));
iblockVector = this_random_pick;
last_4_picks = nan(1,4);
last_4_picks(2) = 6;

while numel(iblockVector)< (n_periodic_trials*length(REGblocks) + n_irreg_trials*length(IRblocks))
    
    % Update scorekeepers
    last_4_picks(2:end) = last_4_picks(1:end-1);
    last_4_picks(1) = this_random_pick;
    
    while this_random_pick==last_4_picks(1)
        if last_4_picks(1) > max(REGblocks) %if last was IRREG, choose REG
            this_random_pick = REGblocks(round(1 + (length(REGblocks)-1).*rand(1,1)));
            
        elseif all(last_4_picks <= max(REGblocks)) %if all 4 last were REG, choose IRREG
            this_random_pick = IRblocks(round(1+(length(IRblocks)-1).*rand(1,1)));
            
        elseif any(last_4_picks <= max(REGblocks)) %if IRREG was one of previous 4 blocks, randomly select which is next
            
            choose_IR = rand(1,1);
            if choose_IR>=0.7 %select IR with 20% chance
                this_random_pick = IRblocks(round(1+(length(IRblocks)-1).*rand(1,1)));
            else
                this_random_pick = REGblocks(round(1 + (length(REGblocks)-1).*rand(1,1)));
            end
            
        end
        
    end
    
    iblockVector = [iblockVector this_random_pick];
    
    % Remove indices of blocks that have reached n trials
    
    
end

findIRs = find(iblockVector==6)


%%

select_periodic_blocks = repmat(REGblocks,1,n_periodic_trials);
select_periodic_blocks = select_periodic_blocks(randperm(length(select_periodic_blocks)));
select_irreg_blocks = repmat(IRblocks,1,n_irreg_trials);
select_irreg_blocks = select_irreg_blocks(randperm(length(select_irreg_blocks)));

% Get number of blocks for each condition
blockNratio = length(select_periodic_blocks)/length(select_irreg_blocks);
[N,D] = rat(blockNratio);

% Interleave regualr and irregular blocks
for isess = 1
    iblockVector = [];
    
    while numel(select_periodic_blocks)>=N && numel(select_irreg_blocks)>=D
        
        % Get indices for irregular blocks
        IRidx=[1,2];
        while any(abs(diff(IRidx))<2) %make sure not 2 consecutive irregular blocks
            IRidx = round(2 + (N+D-2).*rand(1,D));
        end
        
        REGidx = 1:(N+D);
        REGidx(IRidx)=[];
        
        add_blocks = nan(1,N+D);
        add_blocks(REGidx) = select_periodic_blocks(1:N);
        add_blocks(IRidx)  = select_irreg_blocks(1:D);
        
        % Add blocks
        iblockVector = [iblockVector add_blocks];
        
        % Remove
        select_periodic_blocks(1:N) = [];
        select_irreg_blocks(1:D) = [];
        
    end
    
    blockVector = stim_array(iblockVector);
    
end

iblockVector

rateVector



%% Questions

% is the range  2 - 64 Hz  broad enough? (specifically, to hit edges of
% units' rate tuning curves/MTFs)
% Joris fig 9: more or less, yes.


%%


periodic_durations = round(periodic_rates*irreg_chunk_dur).*(1./periodic_rates)

(periodic_durations - irreg_chunk_dur)*1000













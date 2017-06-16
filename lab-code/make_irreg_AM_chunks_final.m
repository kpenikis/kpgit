function make_irreg_AM_chunks_final
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

rng('shuffle')

AMrates = 2.^[1:6];  %2-64 Hz
periodic_rates = [2 4 8 16 32 64];

irreg_chunk_dur = 2*sum(1./AMrates);

try
    
    load IRsequences.mat
    
catch
    
    IRseq_idxs = make_IR_segments;
    
    % Final IR segments
    IR_A = AMrates(IRseq_idxs(1,:));
    IR_B = AMrates(IRseq_idxs(2,:));
    IR_C = AMrates(IRseq_idxs(3,:));
    IR_D = AMrates(IRseq_idxs(4,:));
    
    % save('IRsequences','IR_A','IR_B','IR_C','IR_D','-v7.3')
    
end

IR_rateVec_1     = [IR_A IR_C];
IR_rateVec_1_inv = [IR_C IR_A];
IR_rateVec_2     = [IR_B IR_D];
IR_rateVec_2_inv = [IR_D IR_B];


% Plot the AM rate sequences
figure; hold on
xstart = 0;
plot(xstart+1:xstart+length(IR_rateVec_1),IR_rateVec_1,'+k')
xstart = 10+ xstart+length(IR_rateVec_1);
plot(xstart+1:xstart+length(IR_rateVec_1_inv),IR_rateVec_1_inv,'+b')
xstart = 10 + xstart+length(IR_rateVec_1_inv);
plot(xstart+1:xstart+length(IR_rateVec_2),IR_rateVec_2,'+k')
xstart = 10+ xstart+length(IR_rateVec_2);
plot(xstart+1:xstart+length(IR_rateVec_2_inv),IR_rateVec_2_inv,'+b')
set(gca,'yscale','log')
set(gca,'ytick',AMrates)
xlabel('period number')
ylabel('AM rate (Hz)')



%% Estimating stimulus time, number of trials for a session

ml_target = 2.5;%ml
ml_flowrate = 0.28;%ml/m
target_time = ml_target / ml_flowrate * 60;

n_irreg_trials = 26;
irreg_repeats_time = irreg_chunk_dur * 4 * n_irreg_trials;

n_periodic_rates = numel(periodic_rates);
n_periodic_trials = 26;
periodic_time = (n_periodic_rates+2) * irreg_chunk_dur * n_periodic_trials;

standard_stim_time = irreg_repeats_time + periodic_time
ml_consumed = standard_stim_time/60 * ml_flowrate



%% Piece together long AM vector for continuous sound stream

REGblocks = 1:numel(periodic_rates);
IRblocks  = numel(periodic_rates)+(1:4);
% therefore, blocks are:
%  [ 2 4 8 16 32 64 IR_rateVec_1 IR_rateVec_1_inv IR_rateVec_2 IR_rateVec_2_inv ]
%    1 2 3  4  5  6        7            8                9           10

try
    load Block_Order
    
catch
    REG_pool = repmat(REGblocks([1 2 5 6]),1,n_periodic_trials/2);
    IR_pool  = repmat(IRblocks,1,n_irreg_trials/2);
    
    tic
    
    all_pool = [REG_pool IR_pool];
    SATISFIED = 0;
    
    while SATISFIED==0
        
        all_pool = all_pool(randperm(length(all_pool)));
        
        % conditions:
        % - avoid repeated periodic rates
        % - restrict number of times IR blocks called 3 times in a row
        % - first block must be REG
        consecutiveIRs = find(diff(find(all_pool>=min(IR_pool)),1,2)==1);
        check_REGs = all_pool;
        check_REGs(check_REGs>=min(IR_pool)) = nan;
        if (sum( diff(consecutiveIRs,1,2)==1 ) < 3) && ~any(diff(check_REGs,1,2)==0) && all_pool(1)<min(IR_pool)
            Block_Order = all_pool;
            SATISFIED = 1;
            break
        end
        
    end
    q=toc
    
    Block_Order = [Block_Order Block_Order];
    
    save('Block_Order','Block_Order','-v7.3')
    
end


%% Construct period rate vector


IR_rateVecs = [IR_rateVec_1;...
    IR_rateVec_1_inv;...
    IR_rateVec_2;...
    IR_rateVec_2_inv ];

periodic_nCycles = ceil(periodic_rates*irreg_chunk_dur);
periodic_durations = periodic_nCycles .* (1./periodic_rates);

fullStreamRateVector  = [];
fullStreamBlockVector = [];
track_preREG = [3 4; 4 3; 3 4; 4 3];
for ib = 1:numel(Block_Order)
    
    % REG block
    if Block_Order(ib)<=max(REGblocks)
        
        fullStreamRateVector = [fullStreamRateVector repmat( AMrates(Block_Order(ib)), 1, periodic_nCycles(Block_Order(ib)) ) ];
        fullStreamBlockVector = [fullStreamBlockVector repmat( Block_Order(ib), 1, periodic_nCycles(Block_Order(ib)) ) ];
        
    % IR block
    else
        
        %first add preceding REG block
        thisREG = track_preREG(Block_Order(ib)-6,1);
        fullStreamRateVector = [fullStreamRateVector repmat( AMrates(thisREG), 1, periodic_nCycles(thisREG) ) ];
        fullStreamBlockVector = [fullStreamBlockVector repmat( thisREG, 1, periodic_nCycles(thisREG) ) ];
        track_preREG(Block_Order(ib)-6,:) = fliplr(track_preREG(Block_Order(ib)-6,:));
        
        %then add IR block
        fullStreamRateVector = [fullStreamRateVector IR_rateVecs(Block_Order(ib)-6,:)];
        fullStreamBlockVector = [fullStreamBlockVector repmat( Block_Order(ib), 1, size(IR_rateVecs(Block_Order(ib)-6,:),2) ) ];
        
    end
    
end

figure;
plot(fullStreamRateVector,'-k','LineWidth',1)
hold on
plot(fullStreamBlockVector,'-b','LineWidth',1)
set(gca,'yscale','linear')


%%  Save the vector of rates and a vector of block labels, 

%%%%%%%%%
% Rates
%%%%%%%%%
buffer = [fullStreamRateVector(1) fullStreamRateVector];
% save('fullStream_AMrateVec_25trs_20160323','buffer','-v7.3')

%%%%%%%%%
% Blocks
%%%%%%%%%
stim_array = strsplit(num2str(periodic_rates));
stim_array{end+1} = 'IR_A IR_C';
stim_array{end+1} = 'IR_C IR_A';
stim_array{end+1} = 'IR_B IR_D';
stim_array{end+1} = 'IR_D IR_B';
% 
% buffer=[];
% buffer = [fullStreamBlockVector(1) fullStreamBlockVector];
% save('fullStream_BlockVec_25trs_20160323','buffer','stim_array','-v7.3')


keyboard
aaa=2343;



%% Questions

% is the range  2 - 64 Hz  broad enough? (specifically, to hit edges of
% units' rate tuning curves/MTFs)
% Joris fig 9: more or less, yes.




end










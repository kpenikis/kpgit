function make_IR_abbriged_plusTargetVec
% 

rng('shuffle')

AMrates = 2.^[1:6];  %2-64 Hz
periodic_rates = [2 4 8 16 32 64];

irreg_chunk_dur = 2*sum(1./AMrates);


load IRsequences.mat

IR_rateVec_1     = [IR_A IR_C];
IR_rateVec_1_inv = [IR_C IR_A];
IR_rateVec_2     = [IR_B IR_D];
IR_rateVec_2_inv = [IR_D IR_B];


IR1 = IR_rateVec_1;
IR2 = IR_rateVec_2_inv;
RG1 = 4;
RG2 = 16;
RG3 = 64;

%% Piece together long AM vector for continuous sound stream

REGblocks = [find(periodic_rates==RG1) find(periodic_rates==RG2) find(periodic_rates==RG3)];
IRblocks  = [7 10];
Blocks = [REGblocks IRblocks];
% block IDs are:
%  [ 2 4 8 16 32 64 IR_rateVec_1 IR_rateVec_1_inv IR_rateVec_2 IR_rateVec_2_inv ]
%    1 2 3  4  5  6        7            8                9           10

try 
    load Block_Order_Target
catch
    
    Block_Order=[];
    for ii = 1:40
        
        this_perm = Blocks(randperm(numel(Blocks)));
        while abs(find(this_perm==7) - find(this_perm==10)) ==1
            this_perm = Blocks(randperm(numel(Blocks)));
        end
        
        Block_Order = [Block_Order this_perm];
    end
    
    save('Block_Order_Target','Block_Order','-v7.3')
    
end

%% Construct period rate vector AND target buffer


IR_rateVecs = [IR_rateVec_1;...
    IR_rateVec_2_inv ];
i_IR_switch = [8;8];

periodic_nCycles = ceil(periodic_rates*irreg_chunk_dur);
% periodic_durations = periodic_nCycles .* (1./periodic_rates);

fullStreamRateVector  = [];
fullStreamBlockVector = [];
fullStreamTargetVector = [];

% track_preREG = [3 4; 4 3; 3 4; 4 3];
for ib = 1:numel(Block_Order)
    
    nextTarget = round(rand(1,1));
    if ib>1
        currTarget = fullStreamTargetVector(end);
    else
        currTarget = nextTarget;
    end
    
    
    % REG block
    if Block_Order(ib)<=max(REGblocks)
        
        fullStreamRateVector = [fullStreamRateVector repmat( AMrates(Block_Order(ib)), 1, periodic_nCycles(Block_Order(ib)) ) ];
        fullStreamBlockVector = [fullStreamBlockVector repmat( Block_Order(ib), 1, periodic_nCycles(Block_Order(ib)) ) ];
        
        thisTargetVec = [ repmat( currTarget, 1, floor(0.75*periodic_nCycles(Block_Order(ib))) ) ...
                          repmat( nextTarget, 1, periodic_nCycles(Block_Order(ib)) - floor(0.75*periodic_nCycles(Block_Order(ib))) )];
        fullStreamTargetVector = [fullStreamTargetVector thisTargetVec ];
        
    % IR block
    else
        
        fullStreamRateVector = [fullStreamRateVector IR_rateVecs(Block_Order(ib)==IRblocks,:)];
        fullStreamBlockVector = [fullStreamBlockVector repmat( Block_Order(ib), 1, size(IR_rateVecs(Block_Order(ib)==IRblocks,:),2) ) ];
        
        thisTargetVec = [ repmat( currTarget, 1, i_IR_switch(Block_Order(ib)==IRblocks,1)-1 ) ...
            repmat( nextTarget, 1, size(IR_rateVecs,2) - i_IR_switch(Block_Order(ib)==IRblocks,1)+1 )];
        fullStreamTargetVector = [fullStreamTargetVector thisTargetVec ];
        
    end
    
end

figure; hold on
plot(fullStreamRateVector,'-k','LineWidth',1)
hold on
plot(fullStreamBlockVector,'-b','LineWidth',1)
plot(fullStreamTargetVector,'-r','LineWidth',1)
set(gca,'yscale','linear')


%%  Save the vector of rates and a vector of block labels, 

keyboard

savefolder = 'SwitchSpectrumSession';
if ~exist(savefolder,'dir')
    mkdir(savefolder)
end
addpath(savefolder)

%%%%%%%%%
% Rates
%%%%%%%%%
buffer = [fullStreamRateVector(1) fullStreamRateVector];
save(fullfile(savefolder,'fullStream_AMrateVec_abriged'),'buffer','-v7.3')

%%%%%%%%%
% Blocks
%%%%%%%%%
stim_array = strsplit(num2str(periodic_rates));
stim_array{end+1} = 'IR_A IR_C';
stim_array{end+1} = 'IR_C IR_A';
stim_array{end+1} = 'IR_B IR_D';
stim_array{end+1} = 'IR_D IR_B';

buffer=[];
buffer = [fullStreamBlockVector(1) fullStreamBlockVector];
save(fullfile(savefolder,'fullStream_BlockVec_abriged'),'buffer','stim_array','-v7.3')

%%%%%%%%%%
% Targets
%%%%%%%%%%

buffer=[];
buffer = [fullStreamTargetVector(1) fullStreamTargetVector];
save(fullfile(savefolder,'fullStream_TargetVec_abriged'),'buffer','-v7.3')



aaa=2343;


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




end










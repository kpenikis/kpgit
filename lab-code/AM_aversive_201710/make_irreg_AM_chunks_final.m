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

AMrates = 2.^[1:5];  %2-32 Hz

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


%% Estimating stimulus time, number of trials for a session

irreg_chunk_dur = sum(1./IR_A);

ml_target = 2.5;%ml
ml_flowrate = 0.3;%ml/m
target_time = ml_target / ml_flowrate * 60; %s

n_irreg_trials = 40;
irreg_repeats_time = irreg_chunk_dur * 4 * n_irreg_trials;

n_periodic_rates = numel(AMrates);
n_periodic_trials = 40;
periodic_time = n_periodic_rates * irreg_chunk_dur * n_periodic_trials;

standard_stim_time = irreg_repeats_time + periodic_time
ml_consumed = standard_stim_time/60 * ml_flowrate



%%
%
%  8 IR_A     1
% 16 IR_A     2
%  8 IR_B     3
% 16 IR_B     4
% IR_C 32     5
% IR_C  4     6
% IR_D 32     7
% IR_D  4     8
%  2          9
% 


% Estimate number of trials for each stim transition

running_dur = 0;
Blocks = [];
lastTr = 0; pickTr = 0;
while running_dur < 500 %s
    
    while lastTr == pickTr
        pickTr = ceil(9*rand(1));
    end
    
    Blocks = [Blocks pickTr];
    
    if pickTr<9
        running_dur = running_dur + irreg_chunk_dur+1;
    else
        running_dur = running_dur + 1;
    end
    
    lastTr = pickTr;
    
end

n1 = sum(Blocks==1);
n2 = sum(Blocks==2);
n3 = sum(Blocks==3);
n4 = sum(Blocks==4);
n5 = sum(Blocks==5);
n6 = sum(Blocks==6);
n7 = sum(Blocks==7);
n8 = sum(Blocks==8);
n9 = sum(Blocks==9);



%% Construct period rate vector

REGstim = 1:numel(AMrates);
% therefore, blocks are:
%  [ 2 4 8 16 32 IR_A IR_B IR_C IR_D ]
%    1 2 3  4  5  6    7    8    9  

periodic_nCycles = ceil(AMrates*irreg_chunk_dur);
periodic_durations = periodic_nCycles .* (1./AMrates);


IR_rateVecs = [ IR_A;...
                IR_B;...
                IR_C;...
                IR_D ];
            
StimComboLabels = { '8-IRA';...
                    '16-IRA';...
                    '8-IRB';...
                    '16-IRB';...
                    'IRC-32';...
                    'IRC-4';...
                    'IRD-32';...
                    'IRD-4' };
                    
StimCombos = [ 3 6;...
               4 6;...
               3 7;...
               4 7;...
               8 5;...
               8 2;...
               9 5;...
               9 2 ];          
          
for ist = 1:size(StimCombos,1)
    
    theseStim = StimCombos(ist,:);
    buffer = [];
    
    for iis = 1:2
        if theseStim(iis)<=numel(REGstim) %if periodic
            buffer = [buffer repmat(AMrates(theseStim(iis)), 1,periodic_nCycles(theseStim(iis)) ) ];
        else
            buffer = [buffer IR_rateVecs(theseStim(iis)-numel(REGstim),:)];
        end
        
    end
    buffer = [buffer(1) buffer];
    save(['rateVec_' num2str(ist) '_' StimComboLabels{ist}],'buffer','-v7.3')
    
end

% Lastly save 2 Hz periodic (alone)
buffer = [2 2 2];
save(['rateVec_9_' '2'],'buffer','-v7.3')


end










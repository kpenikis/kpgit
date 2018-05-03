function make_pdcIR_alltransitions
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


%% Remove 64 Hz period from IR blocks

IR_A(IR_A==64) = [];
IR_B(IR_B==64) = [];
IR_C(IR_C==64) = [];
IR_D(IR_D==64) = [];


%% Estimating stimulus time, number of trials for a session

% transtrial_dur = 1 + sum(1./IR_A);
default_dur = 1 + sum(1./IR_A);
double_dur  = 1 + 2*sum(1./IR_A);
equiv_dur   = 1 + 2.5;

iti_dur = 1;
n_reps = 20;

ml_target = 3;%ml
ml_flowrate = 0.28;%ml/m

total_seconds = round( (double_dur+iti_dur) * (4*numel(AMrates) * n_reps) /2);
ml_consumed = total_seconds/60 * ml_flowrate;

fprintf('To get to %i reps per stimulus at a flow rate of %3.2f ml/m, subject must drink for %i seconds and consume %3.3f ml of water.\n',n_reps,ml_flowrate,round(total_seconds),ml_consumed)

target_time = ml_target / ml_flowrate * 60; %s
reps_per_stim = target_time / round( (double_dur+iti_dur) * 2*numel(AMrates) );

fprintf('If subject consumes %3.1f ml of water at a flow rate of %3.2f ml/m, we''ll get ~%i reps per stimulus.\n',ml_target,ml_flowrate,round(reps_per_stim))




keyboard



%%  Jitter rates away from integer ratio relationship

% 
% s = 1;
% JitterVec = ( 2.^ (s.*rand(1,5) - s/2 + [1:5] ) )
% 
% duration = 1000 + round(sum(1000./JitterVec))
% 
% 
% for ir = AMrates
%     IR_A(IR_A==ir) = JitterVec(AMrates==ir);
%     IR_B(IR_B==ir) = JitterVec(AMrates==ir);
%     IR_C(IR_C==ir) = JitterVec(AMrates==ir);
%     IR_D(IR_D==ir) = JitterVec(AMrates==ir);
% end
% 
% 
% IR_E = [16 4 2 32 8];
% IR_F = [2 16 32 8 4];
% IR_G = [32 16 8 4 2];
% IR_H = [2 4 8 16 32];
% for ir = AMrates
%     IR_E(IR_E==ir) = JitterVec(AMrates==ir);
%     IR_F(IR_F==ir) = JitterVec(AMrates==ir);
%     IR_G(IR_G==ir) = JitterVec(AMrates==ir);
%     IR_H(IR_H==ir) = JitterVec(AMrates==ir);
% end
% 
% AMrates = JitterVec;



%% Construct period rate vector

iPdcStim = 1:numel(AMrates);
IR_rateVecs = [ IR_A;...
                IR_B;...
                IR_C;...
                IR_D ];
IR_rateVecLabels = {'A'...
                    'B'...
                    'C'...
                    'D' };

periodic_nCycles = ceil(AMrates*1);
periodic_durations = periodic_nCycles .* (1./AMrates);

default_dur = 1000 + round(sum(1000./IR_A));


% Default IR type

stimfolder = ['AllTransitions_default' num2str(default_dur)];
if ~exist(stimfolder,'dir')
    mkdir(stimfolder);
end
          
for ipd = 1:numel(iPdcStim)
    
    buffer = [];
    
    for iir = 1:size(IR_rateVecs,1)
        
        buffer = [repmat(AMrates(iPdcStim(ipd)), 1, periodic_nCycles(iPdcStim(ipd)) ) IR_rateVecs(iir,:)];
        buffer = [buffer(1) buffer]
                
        save(fullfile(stimfolder,['rateVec_' num2str(round(AMrates(iPdcStim(ipd)))) '_' IR_rateVecLabels{iir}]),'buffer','-v7.3')
        
    end
end


keyboard


% Double IR blocks

IR_E = [16 4 2 32 8];
IR_F = [2 16 32 8 4];
IR_G = [32 16 8 4 2];
IR_H = [2 4 8 16 32];

IR_rateVecs = [ IR_A IR_C;...
                IR_D IR_B;...
                IR_E IR_F;...
                IR_G IR_H ];

IR_rateVecLabels = {'AC'...
                    'DB'...
                    'EF'...
                    'GH' } ;
                
double_dur = 1000+round(sum(1000./IR_rateVecs(1,:)));


stimfolder = ['AllTransitions_double' num2str(double_dur)];
if ~exist(stimfolder,'dir')
    mkdir(stimfolder);
end
          
for ipd = 1:numel(iPdcStim)
    
    buffer = [];
    
    for iir = 1:size(IR_rateVecs,1)
        
        buffer = [repmat(AMrates(iPdcStim(ipd)), 1, periodic_nCycles(iPdcStim(ipd)) ) IR_rateVecs(iir,:)];
        buffer = [buffer(1) buffer]
                
        save(fullfile(stimfolder,['rateVec_' num2str(round(AMrates(iPdcStim(ipd)))) '_' IR_rateVecLabels{iir}]),'buffer','-v7.3')
        
    end
end



% Evivalent time per rate IR blocks

IR_rates = [2 4 4 8 8 8 8 16 16 16 16 16 16 16 16 32 32 32 32 32 32 32 32 32 32 32 32 32 32 32 32];

IR_rateVecLabels = {'W'...
                    'X'...
                    'Y'...
                    'Z' };

equiv_dur = 1000+round(sum(1000./IR_rates));


stimfolder = ['AllTransitions_equivtime_' num2str(equiv_dur)];
if ~exist(stimfolder,'dir')
    mkdir(stimfolder);
end
          
for ipd = 1:numel(iPdcStim)
    
    buffer = [];
    
    for iir = 1:4
        
        buffer = [repmat(AMrates(iPdcStim(ipd)), 1, periodic_nCycles(iPdcStim(ipd)) ) IR_rates(randperm(length(IR_rates)))];
        buffer = [buffer(1) buffer]
                
        save(fullfile(stimfolder,['rateVec_' num2str(AMrates(iPdcStim(ipd))) '_' IR_rateVecLabels{iir}]),'buffer','-v7.3')
        
    end
end




end










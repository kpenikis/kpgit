function make_pdcIR_alltransitions_slow
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


logrates = [1:0.5:4];
AMrates  = 2.^(logrates + 0.2*rand(size(logrates)) -0.1 ) ;  %2-16 Hz

PdDurs   = [481  357  250  166  121  93  67]; %ceil(1000./AMrates)
AMrates  = 1000./PdDurs;

indices = ...
    [...
    4 6 2 1 5 7 3;...
    2 4 7 5 6 1 3];

IR_Y = AMrates(indices(1,:));
IR_Z = AMrates(indices(2,:));

idxPdc  = [1 3 5];
periodic_nCycles = [3 4 8];


%% Estimating stimulus time, number of trials for a session

IR_durations  = [sum(1./IR_Y) sum(1./IR_Z)];
periodic_durations = PdDurs(idxPdc)./1000 .* periodic_nCycles;

iti_dur = 1;
n_reps = 20;

ml_target = 3;%ml
ml_flowrate = 0.26;%ml/m

total_seconds = n_reps * (2*IR_durations(1)+sum(periodic_durations(1:2)) + IR_durations(2)+periodic_durations(3)) + n_reps*IR_durations(2) + 3*n_reps*iti_dur ;
ml_consumed = total_seconds/60 * ml_flowrate;

fprintf('To get to %i reps per stimulus at a flow rate of %3.2f ml/m, subject must drink for %i seconds and consume %3.3f ml of water.\n',n_reps,ml_flowrate,round(total_seconds),ml_consumed)


%% Construct period rate vector

iPdcStim = 1:numel(AMrates);
IR_rateVecs = [ IR_Y;...
                IR_Z ];
IR_rateVecLabels = {'Y'...
                    'Z' };

                
% Default IR type

stimfolder = 'NewTransStimuli';
if ~exist(stimfolder,'dir')
    mkdir(stimfolder);
end
          
for ipd = idxPdc
    
    buffer = [];
    
    switch ipd
        case 1       % 2 to IR_Y 
            buffer = [repmat(AMrates(ipd), 1, periodic_nCycles(ipd==idxPdc) ) IR_rateVecs(1,:)];
            DUR    = ceil(sum(buffer)*1000);
            buffer = [buffer(1) buffer]
            
            save(fullfile(stimfolder,['rateVec_' num2str(round(AMrates(ipd))) '_' IR_rateVecLabels{1} '_' num2str(DUR)]),'buffer','-v7.3')
            
        case 3       % 4 to IR_Y 
            buffer = [repmat(AMrates(ipd), 1, periodic_nCycles(ipd==idxPdc) ) IR_rateVecs(1,:)];
            DUR    = ceil(sum(buffer)*1000);
            buffer = [buffer(1) buffer]
            
            save(fullfile(stimfolder,['rateVec_' num2str(round(AMrates(ipd))) '_' IR_rateVecLabels{1} '_' num2str(DUR)]),'buffer','-v7.3')
            
        case 5       % 8 to IR_Z 
            buffer = [repmat(AMrates(ipd), 1, periodic_nCycles(ipd==idxPdc) ) IR_rateVecs(2,:)];
            DUR    = ceil(sum(buffer)*1000);
            buffer = [buffer(1) buffer]
            
            save(fullfile(stimfolder,['rateVec_' num2str(round(AMrates(ipd))) '_' IR_rateVecLabels{2} '_' num2str(DUR)]),'buffer','-v7.3')
            
    end
    
end

% IR_Z alone
buffer = [];
buffer = IR_rateVecs(2,:);
DUR    = ceil(sum(buffer)*1000);
buffer = [buffer(1) buffer]
save(fullfile(stimfolder,['rateVec_' IR_rateVecLabels{2} '_' num2str(DUR)]),'buffer','-v7.3')

% 4 hz alone for ITI (and Warn)
buffer = [];
buffer = repmat(AMrates(ipd), 1, periodic_nCycles(ipd==idxPdc) ) ;
DUR    = ceil(sum(buffer)*1000);
buffer = [buffer(1) buffer]

save(fullfile(stimfolder,['rateVec_' num2str(round(AMrates(ipd))) '_' num2str(DUR)]),'buffer','-v7.3')




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
                
IR_durations = 1000+round(sum(1000./IR_rateVecs(1,:)));


stimfolder = ['AllTransitions_double' num2str(IR_durations)];
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


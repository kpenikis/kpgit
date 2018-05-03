function convert_stream_to_trials

%  [ 2 4 8 16 32 64 IR_rateVec_1 IR_rateVec_1_inv IR_rateVec_2 IR_rateVec_2_inv ]
%    1 2 3  4  5  6        7            8                9           10

% IR_rateVec_1     = [IR_A IR_C];
% IR_rateVec_1_inv = [IR_C IR_A];
% IR_rateVec_2     = [IR_B IR_D];
% IR_rateVec_2_inv = [IR_D IR_B];


savefolder = 'REG-IR_trials';
if ~exist(savefolder,'dir')
    mkdir(savefolder)
end
addpath(savefolder)

rateBuffer  = load('fullStream_AMrateVec_25trs_20160323.mat');
rateBuffer  = rateBuffer.buffer;
blockBuffer = load('fullStream_BlockVec_25trs_20160323.mat');
blockBuffer = blockBuffer.buffer;

transition_indices = find(diff(blockBuffer)); %last period of previous block
transition_indices = [0 transition_indices];

for ib = 2:numel(transition_indices)
    
    this_block  = blockBuffer(transition_indices(ib-1)+1:transition_indices(ib));
    if ~all(diff(this_block)==0)
        keyboard
    end
    if this_block(1)<=6
        blktype = 'REG';
    elseif this_block(1)>=7
        blktype = 'IR';
    end
    
    buffer = rateBuffer(transition_indices(ib-1)+1:transition_indices(ib));
%     buffer = [these_rates(1) these_rates];
    
%     eval(sprintf('blk%i_%s_%i = buffer;',(ib-1),blktype,this_block(1)))
    
%     save(fullfile(savefolder,sprintf('blk%i_%s_%i',(ib-1),blktype,this_block(1))),'-v7.3')
    
    clear buffer 

end


end


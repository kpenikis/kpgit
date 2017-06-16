function set_Target_indices

%  [ 2 4 8 16 32 64 IR_rateVec_1 IR_rateVec_1_inv IR_rateVec_2 IR_rateVec_2_inv ]
%    1 2 3  4  5  6        7            8                9           10

load('IRsequences.mat')
IR_rateVec_1     = [IR_A IR_C];
IR_rateVec_1_inv = [IR_C IR_A];
IR_rateVec_2     = [IR_B IR_D];
IR_rateVec_2_inv = [IR_D IR_B];

rateBuffer  = load('fullStream_AMrateVec_25trs_20160323.mat');
rateBuffer  = rateBuffer.buffer;
blockBuffer = load('fullStream_BlockVec_25trs_20160323.mat');
blockBuffer = blockBuffer.buffer;

transition_indices = find(diff(blockBuffer)); %last period of previous block
transition_indices = [0 transition_indices];

%%


TargetTriggers = zeros(size(rateBuffer));
trans_lastblock = 1;
hold2  = false;
hold5  = false;
hold8  = false;
hold10 = false;

for ib = 2:numel(transition_indices)
    
    this_block  = blockBuffer(transition_indices(ib-1)+1:transition_indices(ib));
    if ~all(diff(this_block)==0)
        keyboard
    end
        
    if (this_block(1)==8 && ~hold8) || (this_block(1)==10 && ~hold10)
        
        if rand(1,1)>0 && trans_lastblock == 0
            TargetTriggers( transition_indices(ib-1)+1 + 0.75*(numel(this_block)) ) = 1;
            trans_lastblock = 1;
        else
            trans_lastblock = 0;
        end
        
        
    elseif ((this_block(1)==2 && ~hold2) || (this_block(1)==5 && ~hold5)) && rand(1,1)>0.1
        
        if rand(1,1)>0 && trans_lastblock == 0
            TargetTriggers( transition_indices(ib-1)+1 + ceil(0.75*(numel(this_block))) ) = 1;
            trans_lastblock = 1;
        else
            trans_lastblock = 0;
        end
    end
    
    hold2  = sum(blockBuffer(TargetTriggers==1)==2)  > 11;
    hold5  = sum(blockBuffer(TargetTriggers==1)==5)  > 11;
    hold8  = sum(blockBuffer(TargetTriggers==1)==8)  > 11;
    hold10 = sum(blockBuffer(TargetTriggers==1)==10) > 11;
    
end

sum(TargetTriggers)
sum(blockBuffer(TargetTriggers==1)==2)
sum(blockBuffer(TargetTriggers==1)==5)
sum(blockBuffer(TargetTriggers==1)==8)
sum(blockBuffer(TargetTriggers==1)==10)



%%
% 
% buffer = TargetTriggers;
% save('fullStream_TargetIndices_2ea_20160323','buffer','-v7.3')
% 
% 
% blk_indices = find(eachblock==4);
% eachblock(blk_indices+1)


end



function create_AM_jitter_matfiles_uniform_log2_v3
% New version for creating AM jitter stimuli. 
% Now done formulaically, along presumed perceptual axis.
% Also now always includes the "centroid" AM rate in the vectors, currently
% designed to be in the center position (i.e. 4th rate of vector length 7).
% This allows a standard comparison across randomized and ranked stimuli.
% 
% KP, 2016.

rng('shuffle');


%Stimulus parameters
AMrate      = 8; %Hz
AMrate_exp = log2(AMrate);
log2steps   = [0 0.1 0.2 0.3 0.4 0.5 0.7 1 1.5 2] ; 

N_periods   = 7; %must be an odd number, because "centroid" rate will always be included
AMrate_pos  = (N_periods+1)/2;
uniqueRates = N_periods;
reps = 5;


%Saving info
pn = 'D:\stim\AMjitter';
savedir = sprintf('AM_%iHz_final',AMrate);

%Plotting info
p_cols = hsv(length(log2steps));
figure(1); clf
plot([1 uniqueRates+1],[AMrate AMrate],'-k')

for ii = 1:numel(log2steps)
    
    rateMin = 2^(AMrate_exp - log2steps(ii) );
    rateMax = 2^(AMrate_exp + log2steps(ii) );
        
    %Get rate values equally spaced in logarithmic axis
    if log2steps(ii)>0
        
        steps = (uniqueRates+1)/2;
        
        rateVec = [2.^ linspace( AMrate_exp-log2steps(ii), AMrate_exp , steps)...
                   2.^ linspace( AMrate_exp, AMrate_exp+log2steps(ii), steps)];
               
        rateVec(steps:steps+1)=[]; %remove repeated rates. they will be added later at a fixed position.
        
        
        %Now permute rates for randomized vectors
        rVidx=[]; buffer=[];
        
        if N_periods > uniqueRates      % long vectors for Appetitive training
            
            rateVec = [rateVec AMrate]; %add centroid rate, order doesn't matter
            
            while numel(rVidx)<N_periods
                rVidx = [rVidx randperm(length(rateVec))];
            end
            
            buffer = rateVec(rVidx(1:N_periods));
            
            filename = sprintf('AM_%iHz_%s_rnd.mat',AMrate,num2str(log2steps(ii)*100));
            savefilename = fullfile(pn,'AppetitiveTraining',filename);
            save(savefilename,'buffer','-v7.3')
            
            
        elseif N_periods == uniqueRates  % short vectors for ephys and behavior
            
            %Save ranked versions first
            
             %increasing rate
            buffer = sort([rateVec AMrate],'ascend');
            buffer = [buffer(1) buffer];
            
            filename = sprintf('AM_%iHz_%s_inc.mat',AMrate,num2str(log2steps(ii)*100));
            savefilename = fullfile(pn,savedir,filename);
%             save(savefilename,'buffer','-v7.3')
            
             %decreasing rate
            buffer = sort([rateVec AMrate],'descend');
            buffer = [buffer(1) buffer];
            
            filename = sprintf('AM_%iHz_%s_dec.mat',AMrate,num2str(log2steps(ii)*100));
            savefilename = fullfile(pn,savedir,filename);
%             save(savefilename,'buffer','-v7.3')
            
            
            %Shuffle indices to get random rate vectors
            for rep = 1:reps
                
                buffer=[];
                rVidx = randperm(length(rateVec));
                
                buffer = [rateVec(rVidx(1)),...   %first rate will be skipped in program, but is used for TTL delay
                          rateVec(rVidx(1:AMrate_pos-1)),...
                          AMrate,...   %insert centroid AM rate in desired (middle) position
                          rateVec(rVidx(AMrate_pos:end))];
                
                
                %Save randomized rate vector file for each trial
                
                filename = sprintf('AM_%iHz_%s_rnd-%i.mat',AMrate,num2str(log2steps(ii)*100),rep);
                savefilename = fullfile(pn,savedir,filename);
%                 save(savefilename,'buffer','-v7.3')
                
                %Plot rate vectors
%                 if log2steps(ii)==1 
                    hold on
                    plot( 1:N_periods+1 , buffer , 'o-', 'Color',p_cols(ii,:) )
                    plot([1 1; N_periods+1 N_periods+1],...
                        [rateMin rateMax; rateMin rateMax],...
                        '--','Color',p_cols(ii,:))
                    hold off
%                 end
                
            end
        end
        
    elseif log2steps(ii)==0 %create and save vector for strictly periodic case
        
        buffer = repmat(AMrate,[1 N_periods+1]);
        
        filename = sprintf('AM_%iHz_%s.mat',AMrate,num2str(log2steps(ii)*100));
        if N_periods > uniqueRates
            savefilename = fullfile(pn,'AppetitiveTraining',filename);
        elseif N_periods == uniqueRates
            savefilename = fullfile(pn,savedir,filename);
        end
%         save(savefilename,'buffer','-v7.3')
        
    end
    
end
%Finish plot
set(gca,'YScale','log')
ylim([2.^(AMrate_exp-log2steps(ii)-0.2) 2.^(AMrate_exp+log2steps(ii)+0.2)])

end








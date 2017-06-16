function [x,y,prev250FR,prev100FR,prev50FR,totalFR] = standardPd_getspikes(data)
% 

global fn subject

% Get vector of rates for this stimulus
stimdir = fullfile(fn.raw,subject,sprintf('Block-%i_Stim',data.block(1)));
rateVec = load(fullfile(stimdir,data.stimfn));
rateVec = rateVec.buffer;
rateVec = rateVec(2:end);

ntr = length(data.tr_idx);

pd_standard = ceil((length(rateVec))/2);
if pd_standard ~=4, keyboard, end

% Get time window of the standard period
t0 = ceil(data.AMonset + 0.75*1000/rateVec(1) + sum(1000./rateVec(2:pd_standard-1)));
tVec = [t0 t0+1000/rateVec(pd_standard)-1];

% Get spiketimes
x = data.x;
y = data.y;

% Get FR in time window 250 & 100 ms prior to standard period
prev250FR = sum( data.x>=(t0-250) & data.x<=t0 )/ntr/250*1000;
prev100FR = sum( data.x>=(t0-100) & data.x<=t0 )/ntr/100*1000;
prev50FR  = sum( data.x>=(t0-50)  & data.x<=t0 )/ntr/50*1000;
totalFR   = sum( data.x>=0        & data.x<=t0 )/ntr/t0*1000;

% Restrict spikes to 4 Hz standard period window only
y ( data.x<tVec(1) | data.x>tVec(2) ) = [];
x ( data.x<tVec(1) | data.x>tVec(2) ) = [];

% Correct spiketimes to timestamp of period start
x = x - tVec(1) + 1;

end

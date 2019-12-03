function levelSignal = sound2dB_AM(ampSignal,dBSPL)
% levelSignal = sound2dB_AM(ampSignal,dBSPL)
%  
%  In AM sessions, sound signal was saved before calibration macro. Thus,
%  first adjust as done in the macro, then convert to dB scale. 
%
%  KP, 2019-12

levelSignal = log10( 1.7402 * 10^((dBSPL-94)/20) * ampSignal ) * 20;
levelSignal = real(levelSignal);

end

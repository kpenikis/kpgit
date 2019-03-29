% epDataPhys_to_DAT
%  Takes the phys data from epData and creates a binary (.dat) file of the 
%  raw ephys data, the format required for KiloSort. 
%  KP, 2017-09

fidn=fopen(fullfile(saveDataDir,datName),'w');
fwrite(fidn,int16(epData.streams.Wave.data*1e6),'int16'); %(int16 format requires integers)
fclose(fidn);

function write_Phys_to_DAT(fwpath,Phys,fwname)
% write_Phys_to_DAT(filewritepath,Phys)
%  Takes the Phys data structure and creates a binary (.dat) file of the 
%  raw ephys data, the format required for KiloSort. 
%    fwpath: folder to save the dat file in
%    Phys:   raw data matrix (channels x samples)
%      **TO DO: check running with unfiltered data straight from epData
%    fwname: name for new .dat file (by default, will be the same as the
%    Phys filename for now)
%  Option to run manually, with no input arguments, and will be prompted to
%  select the write path and the Phys data structure in a UI. 
%
%  KP, 2017-09.
%

% Select Phys/raw data file if not given as an input
if nargin<2
    fn = set_paths_directories;
    [filename,pathname] = uigetfile('*Phys.mat','Select a Phys file',fn.processed);
    load([pathname filename])
    fwname = filename(1:end-4);
end
% Select write path if not given as an input
if nargin<1
    fwpath = uigetdir('/Users/kpenikis/Documents/SanesLab/LocalData');
end

% Write the .dat file
fidn=fopen(fullfile(fwpath,[fwname '.dat']),'w');
fwrite(fidn,int16(Phys*1e6),'int16'); %(int16 format requires integers)
fclose(fidn);

end
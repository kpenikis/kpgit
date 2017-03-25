function Data = downsize_Data_structure(Data)

addpath('/Users/kpenikis/Documents/MATLAB/psignifit_git')


subject = 'IIIf_230115';
session = 'LA';
datadir  = '/Users/kpenikis/Documents/SanesLab/Data/AMJitter/ProcessedData';
filename = sprintf( '%s_sess-%s_Data',subject,session); 
if nargin<1
    disp('loading Data file')
    load(fullfile(datadir,subject,filename));
    disp('reformatting...')
end

% Data_test = struct;

for ich = 1:numel(Data.ch)
    
    clear d
    
    for iclu = 1:numel(Data.ch(ich).clu)
        
        d.clu(iclu).label  = Data.ch(ich).clu(iclu).label;
        d.clu(iclu).raster = Data.ch(ich).clu(iclu).raster;
        d.clu(iclu).stim   = Data.ch(ich).clu(iclu).analyses.stim;
        Data.ch(ich).clu(iclu).stim = Data.ch(ich).clu(iclu).analyses.stim;
        
        for is = 1:numel(Data.ch(ich).clu(iclu).stim)
            for ip = 1:numel(Data.ch(ich).clu(iclu).stim(is).pars)
                
                for iVt = {'nm_jitter' 'nm_depth'}
                    indVar_type = iVt{:};
                    if isfield(Data.ch(ich).clu(iclu).stim(is).pars(ip),indVar_type)
                        
                        classifier = struct;
                        for icp = 1:numel(Data.ch(ich).clu(iclu).stim(is).pars(ip).(indVar_type).classifier)
                            
                            classifier(icp).metric  = Data.ch(ich).clu(iclu).stim(is).pars(ip).(indVar_type).classifier(icp).params.metric;
                            classifier(icp).binsize = Data.ch(ich).clu(iclu).stim(is).pars(ip).(indVar_type).classifier(icp).params.binsize;
                            
                            clear output
                            output = Data.ch(ich).clu(iclu).stim(is).pars(ip).(indVar_type).classifier(icp).output;
                            fns = fieldnames(output);
                            for ifn = 1:numel(fns)
                                
                                if any(strcmp(fns{ifn},{'PYdata' 'dprime'}))
                                    for icond = 1:numel(output.(fns{ifn}))
                                        clear input
                                        try
                                        input.(fns{ifn})(:,:,icond) = output.(fns{ifn}){icond};
                                        catch
                                            keyboard
                                        end
                                    end
                                    classifier(icp).(fns{ifn}) = input.(fns{ifn});
                                    
                                elseif strcmp(fns{ifn},'stim') && icp==1 && strcmp(indVar_type,'nm_jitter')
                                    d.clu(iclu).stim(is).pars(ip).stim_jitter = output.(fns{ifn});
                                    output = rmfield(output,fns{ifn});
                                elseif strcmp(fns{ifn},'stim') && icp==1 && strcmp(indVar_type,'nm_depth')
                                    d.clu(iclu).stim(is).pars(ip).stim_depth = output.(fns{ifn});
                                    output = rmfield(output,fns{ifn});
                                
                                elseif strcmp(fns{ifn},'fitdata')
                                    %remove some extra fields from psignifit output
                                    try
                                        for ipy = 1:numel(output.(fns{ifn}).PYdata)
                                            output.(fns{ifn}).PYdata(ipy).results = rmfield(output.(fns{ifn}).PYdata(ipy).results,{'X1D' 'marginals' 'marginalsX' 'marginalsW'});
                                        end
                                    catch
                                        keyboard
                                    end
                                    classifier(icp).(fns{ifn}) = output.(fns{ifn});
                                end
                                
                            end
                        end
%                         d.(varname).clu(iclu).stim(is).pars(ip).(indVar_type) = rmfield(Data.ch(ich).clu(iclu).stim(is).pars(ip).(indVar_type),'classifier');
                        d.clu(iclu).stim(is).pars(ip).(indVar_type).classifier = classifier;
                        
                        d.clu(iclu).stim(is).pars(ip).(indVar_type).formulaFR = Data.ch(ich).clu(iclu).stim(is).pars(ip).(indVar_type).formulaFR.output;
                        d.clu(iclu).stim(is).pars(ip).(indVar_type).formulaFR = rmfield(d.clu(iclu).stim(is).pars(ip).(indVar_type).formulaFR,'stim');
                        
                    end
                end

            end %ip
        end %is
        
    end %iclu
    
    if exist('d','var')
        eval(sprintf('ch%i = d;',ich));
    end
    
end


filename = sprintf( '%s_sess-%s_Data',subject,session); 
disp('saving...')
save(fullfile(datadir,subject,filename),'-regexp','^ch','-v7.3');
save(fullfile(datadir,subject,filename),'Data','-v7.3');
disp('done.')



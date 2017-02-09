function calc_metrics(directoryname,figuredirectory,n,animalIDs)
%calc_metrics(directoryname,n)
%
%Calculates firing rate, power, associated d's and d-based thresholds for
%each cluster. Thresholds are determined by fitting dprime neurometric
%functions with a sigmoid and finding the depth at which the fit == 1.
%Site is considered sound responsive if at least one of the fits returned a
%valid threshold.
%
%Input variables:
%   directoryname:directory containing data to analyze, organized by day,
%                then by animal
%
%   n: seconds after stimulus onset over which data should be analyzed
%
%   animalIDs: cell array of animal IDs to be included in the analysis
%
%
%Written by ML Caras Dec 5 2015

warning('off','stats:nlinfit:ModelConstantWRTParam')
warning('off','stats:nlinfit:IllConditionedJacobian')
warning('off','MATLAB:rankDeficientMatrix')

set(0,'DefaultTextInterpreter','none');
fs = 500; %samples for FFT analysis

%Find list of day folders
[days,dayIndex]= findRealDirs(directoryname);
days = days(dayIndex);

for which_day = 1:numel(days)
    
    %Select one day and send that day's data to a single instance on the HPC
    dayname = [directoryname,days(which_day).name,'/'];
    
    %Get a list of folders in the directory (each folder == one type of data)
    [folders,folderIndex]= findRealDirs(dayname);
    folders = folders(folderIndex);
    
    ind = find_index(folders,'combined');
    
    combined_folder = [dayname,folders(ind).name];
    
    %Get a list of folders in the directory (each folder == one animal)
    [animals,animalIndex] = findRealDirs(combined_folder);
    animals = animals(animalIndex);
    
    %For each animal, load the data and do some calculations.
    for which_animal = 1:numel(animals)
        
        %Define animal subfolder
        animal_folder = [combined_folder,'/',animals(which_animal).name];
        ID = animals(which_animal).name(end-5:end);
        
        %Make sure animal is one of the ones we want to analyze
        if ~any(ismember(animalIDs,ID))
            continue
        end
        
        %Get list of files in the directory
        [files,fileIndex] = listFiles(animal_folder,'*.mat');
        files = files(fileIndex);
        
        %For each file...
        for which_file = 1:numel(files)
            filename = files(which_file).name;
            data_file = [animal_folder,'/',filename];
            
            %Load file
            load(data_file);
            
            %Pull out file parameters
            MF = double(data.behavior.trial_log.fm(1)); %Modulation frequency (Hz)
            duration = data.behavior.trial_log.trial_duration(1);%Duration of stim (sec)
            
            %Pull out stimulus values and the number of trials each
            %stimulus was presented
            [numtrials,stims] = grpstats(data.behavior.PROCESSED.stim_log,...
                data.behavior.PROCESSED.stim_log,{'numel','gname'});
            stims = cell2mat(cellfun(@str2num,stims,'uniformoutput',0));
            stims = make_stim_log(stims);
            
            %Pull out channels
            channels = fieldnames(data.physiology);
            

            %For each channel...
            for which_channel = 1:numel(channels)
                
                %Pull out channel data
                channeldata = data.physiology.(channels{which_channel});
                
                %For each cluster...
                for which_cluster = 1:numel(channeldata.clusters)
                    
                    %Pull out cluster data
                    cluster_data = channeldata.clusters(which_cluster);
                    
                    %Skip cluster if it's not a multi unit or a single unit
                    if ~any(strcmp(cluster_data.clusterType{1},{'multi-unit','good unit'}))
                        continue
                    end
                    
                    %Initialize empty matrices and figure
                    FRmat = [];
                    Powermat = [];
                    f1 = myplot;
                                      
                    %Restrict spikes to those that occurred within n
                    %seconds of stimulus onset
                    %Datamatrix = [trial, stimulus, spiketime (s)]
                    cluster_data.datamatrix = cluster_data.datamatrix(cluster_data.datamatrix(:,3) <=n,:);
                    
                    %For each stimulus
                    for which_stim = 1:numel(stims)
                        n_trials = numtrials(which_stim);
                        stimdata = cluster_data.datamatrix(cluster_data.datamatrix(:,2) == stims(which_stim),:);
                        
                        %--------------------------------
                        %CALCULATE FIRING RATE
                        %---------------------------------
                        FRmat = calculateFR(stimdata,n,n_trials,FRmat,stims,which_stim);
                        
                        %--------------------------------
                        %CALCULATE POWER
                        %--------------------------------
                        Powermat = calculatepower(Powermat,stimdata,...
                            n_trials,stims,which_stim,MF,fs,duration);
                        
                    end
                    
                    %------------------------------------
                    %CALCULATE FR-BASED d' VALUES
                    %------------------------------------
                    dprime_mat = calculatedprime(FRmat);
                    
                    %------------------------------------
                    %CALCULATE POWER-BASED d' VALUES
                    %------------------------------------
                    Power_dprime_mat = calculatedprime(Powermat);
                    
                    
                    %------------------------------------
                    %FIND FR dprime-BASED THRESHOLD
                    %------------------------------------
                    [xFRfit,yFRfit,dprime_threshold,dprime_threshold_15] = ...
                        calculate_dprime_threshold(dprime_mat,numtrials,stims,1);
                    
                    %------------------------------------
                    %FIND Power dprime-BASED THRESHOLD
                    %------------------------------------
                    [xPowerfit,yPowerfit,Power_dprime_threshold,Power_dprime_threshold_15] = ...
                        calculate_dprime_threshold(Power_dprime_mat,numtrials,stims,2);
                    
                    
                    
                    %If we plotted data
                    if ~isempty(get(f1,'children'))
                        
                        %Add a title
                        depth = num2str(data.behavior.metadata.recording_depth);
                        date = data.behavior.metadata.session_date;
                        id = data.behavior.metadata.nyu_id;
                        ch = channels{which_channel};
                        clusType = cluster_data.clusterType{1};
                        clusID = num2str(cluster_data.clusterID);
                        
                        
                        suptitle([id,' ', date, ' ' ,depth, 'um ',ch,' ',clusType, ' cluster ',clusID],gcf)
                        
                        
                        %Save the figure
                        depth = num2str(str2num(depth)*10);
                        figname = [filename(1:end-4),'_',depth,'_',ch,'_',num2str(which_cluster)];
                        set(f1,'PaperPositionMode','auto');
                        print(f1,'-painter','-depsc', [figuredirectory,figname])
                        
                    end
                    
                    close all
                    
           
                    %Append FR mat to cluster
                    data.physiology.(channels{which_channel}).clusters(which_cluster).FRmat = FRmat;
                    
                    %Append Power mat to cluster
                    data.physiology.(channels{which_channel}).clusters(which_cluster).Powermat = Powermat;
                    
                    %Append FRdprime mat to cluster
                    data.physiology.(channels{which_channel}).clusters(which_cluster).FR_dprime_mat = dprime_mat;
                    
                    %Append Powerdprime mat to cluster
                    data.physiology.(channels{which_channel}).clusters(which_cluster).Power_dprime_mat = Power_dprime_mat;
                    
                    %Append FRdprime threshold to cluster (d' == 1)
                    data.physiology.(channels{which_channel}).clusters(which_cluster).FR_dprime_threshold = dprime_threshold;
                    
                    %Append Powerdprime threshold to cluster (d' == 1)
                    data.physiology.(channels{which_channel}).clusters(which_cluster).Power_dprime_threshold = Power_dprime_threshold;
                    
                    %Append FRdprime threshold to cluster (d' == 1.5)
                    data.physiology.(channels{which_channel}).clusters(which_cluster).FR_dprime_threshold_15 = dprime_threshold_15;
                    
                    %Append Powerdprime threshold to cluster (d' == 1.5)
                    data.physiology.(channels{which_channel}).clusters(which_cluster).Power_dprime_threshold_15 = Power_dprime_threshold_15;

                    %Append FRdprime fit to cluster
                    data.physiology.(channels{which_channel}).clusters(which_cluster).FR_dprime_fit = [xFRfit',yFRfit'];
                    
                    %Append Powerdprime fit to cluster
                    data.physiology.(channels{which_channel}).clusters(which_cluster).Power_dprime_fit = [xPowerfit',yPowerfit'];
                    
                    %Append timescale of analysis to cluster
                    data.physiology.(channels{which_channel}).clusters(which_cluster).timescale = n;
                    
                end

            end
            
            
            %Save data file
            save(data_file,'data','-append')
            disp(['Saved ' data_file])

            
        end
        
        
        
    end
    
    
    
end



end


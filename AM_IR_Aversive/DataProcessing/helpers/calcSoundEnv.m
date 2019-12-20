function ENV = calcSoundEnv(InputData,fs,ENV_type)
% calcSoundEnv(InputData,ENV_type)
% 
%  Methods from extracting speech envelope for vocoding. Now use for post
%  hoc analyses too. 
% 
%  KP, 2019 


%%

data = InputData;%/max(InputData);
data = abs(data) - min(abs(data));

if nargin<3
    ENV_type = 'HilbLP';
end

minpeakdist_ms    = 12;
minpeakdist_samps = round(minpeakdist_ms/1000 * fs);


switch ENV_type
    
    case 'Peaks'
        %~~~~~~~~~~~~~~~~~~~
        % ENVELOPE METHOD 1: find peaks, must be separated by at least X samples
        
        
        ENV = envelope( data + (rand(size(data))*0.006-0.003) ,minpeakdist_samps,'peak');
        ENV(ENV<0) = 0;
        
%         ENV   = ENV/max(ENV);
        
        
    case 'HilbLP'
        %~~~~~~~~~~~~~~~~~~~
        % ENVELOPE METHOD 2: Hilbert analytic signal, then low pass filtered
        
%         hilbtrans = abs(hilbert(data));
        hilbtrans = real(hilbert(data));
        
        % Lowpass
%         lpFilt = designfilt('lowpassfir','PassbandFrequency',45, ...
%             'StopbandFrequency',60,'DesignMethod','kaiserwin','SampleRate',fs);
        lpFilt = designfilt('lowpassiir','FilterOrder',2, ...
            'PassbandFrequency',45,'PassbandRipple',0.2,'SampleRate',fs);
%         fvtool(lpFilt)
        
        ENV = filtfilt(lpFilt,hilbtrans);
        ENV = filtfilt(lpFilt,ENV);
        ENV(ENV<0) = 0;
%         ENV = ENV/max(ENV);
        
        
    case 'RMS'
        %~~~~~~~~~~~~~~~~~~~
        % ENVELOPE METHOD 3: rms
        
        ENV = envelope(data,minpeakdist_samps,'rms');
        ENV(ENV<0) = 0;
%         ENV = ENV/max(ENV);
        
end


% %% Plot data and envelopes
% 
% L = length(data);
% TimeVec = linspace(0,L/fs,L);
% 
% % Plot
% figure;
% plot(TimeVec,data,'k')
% hold on
% plot(TimeVec,ENV,'LineWidth',1.5)
% 
% xlim([0 max(TimeVec)])
% ylim([-1 1.3])


end

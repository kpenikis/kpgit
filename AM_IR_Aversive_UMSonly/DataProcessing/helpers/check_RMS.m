function isamp = check_RMS( data, fs )
% Remove beginning of session if noisy, ie if started before headstage 
% turned on or animal was put in room. Should help sorting, because
% threshold gets messed up if a large portion of the session is noisier. 
% 
% KP, 2018-12.


fprintf(' removing noisy pre-session portion...')

begRMS = rms(data(1,1:round(fs*10))); % first 10 s
endRMS = rms(data(1,(end-round(fs*21)):(end-round(fs*1)) )); % last 20 s

if (begRMS-1.2*endRMS)<0.01 %while testing
    keyboard
end

sampwin  = round(fs*2); %2000 ms
sampstep = round(fs/10); %100 ms

isamp = 1;
RMS_HIGH = 1;
while RMS_HIGH
    thisRMS = rms(data(1,isamp:(isamp+sampwin-1)));
    if thisRMS < 1.2*endRMS
        RMS_HIGH = 0;
    else
        isamp = isamp+sampstep;
    end
end

% Plot to check
begwins  = 1:sampstep:size(data,2);
roughRMS = nan(size(begwins));
for iw = 1:numel(begwins)
    roughRMS(iw) = rms(data(1,begwins(iw):min(begwins(iw)+sampwin-1,size(data,2))));
end

hf=figure;
plot(begwins,roughRMS,'k')
hold on
plot([isamp isamp],[0 begRMS],'r','LineWidth',2)

keyboard
close(hf);

end

function envelope_to_Wav(ENV,car_type,fs,savepathname)

% Make sure envelope is in range [0 1]
Envelope = ENV - min(ENV);
Envelope = Envelope/max(Envelope);


% Create carrier signal
switch car_type
    case 'WN'
        Carrier = rand(size(ENV))*2-1;
end

data = Envelope .* Carrier;
data = data/max(data);

audiowrite(savepathname,data,fs)

end
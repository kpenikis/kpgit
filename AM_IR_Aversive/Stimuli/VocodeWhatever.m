
addpath('Stimuli')

modul = audioread('/Volumes/GoogleDrive/My Drive/Sanes/DATADIR/AMaversive/Stimuli/SpeechRecordings/VocodedSpeechStimuli/pigeons are sent by the government 3_01.wav');
carrier = audioread('/Volumes/GoogleDrive/My Drive/Sanes/DATADIR/AMaversive/Stimuli/SpeechRecordings/VocodedSpeechStimuli/Hallelujah - audio_01.wav');
carrier = rand(size(modul));

chan    = 1;
numband = chan;
overlap = 1/10;

y = WiscChanVocoder(carrier, modul, chan, numband, overlap);

figure;
% plot(envelope(y,10,'rms'))
plot(y)

sound(y)

audiowrite('/Volumes/GoogleDrive/My Drive/Sanes/DATADIR/AMaversive/Stimuli/SpeechRecordings/VocodedSpeechStimuli/pigeons vocoded noise 1ch.wav',y,8192);


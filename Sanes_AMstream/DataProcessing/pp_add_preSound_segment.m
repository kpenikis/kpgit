function SoundData = pp_add_preSound_segment(SoundData,Info,Phys)

ds = round(Info.fs/Info.fs_sound);

% find onset of unmodulated noise
unmodOnset = find(diff(SoundData(8,:))==11);

figure; plot(Phys(8, ds.*[1 : unmodOnset(1)]))

% plot
figure;
plot( Phys(8, ds.*[1 : unmodOnset(1)+3*round(Info.fs_sound)])./max(Phys(1, ds.*[1 : unmodOnset(1)+3*round(Info.fs_sound)])), 'Color',[0.5 0.5 0.5])
xlim([1 unmodOnset(1)+3*round(Info.fs_sound)])
hold on
plot( SoundData(8, 1:(unmodOnset(1)+3*round(Info.fs_sound)))./max(SoundData(8, 1:(unmodOnset(1)+3*round(Info.fs_sound)))), 'LineWidth',3 )











end







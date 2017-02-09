
seconds = 1:1200;

MLperMin_Steady = 0.3;

consumed_Steady = MLperMin_Steady/60*seconds;

figure(1); clf
plot(seconds,consumed_Steady)


flowRate_Pulse = 0.6;
timeOn  = 0.3;
timeOff = 1-timeOn;

MLperMin_Pulse = timeOn*flowRate_Pulse;
consumed_Pulse = MLperMin_Pulse/60*seconds;

hold on
plot(seconds,consumed_Pulse)

xlabel('Time in passive experiment')
ylabel('ML water consumed')





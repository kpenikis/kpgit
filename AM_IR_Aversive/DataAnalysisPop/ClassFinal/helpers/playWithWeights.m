beta=([   -0.1025
-0.0686
0.1214
0.5951
-0.2087
0.1629
-0.1841
-0.0667
-0.8035
-0.7429
-0.4563
-0.2279
-0.2216
-0.6020
-0.5519
-0.4204
-0.9343
-0.3512
0.0329
-0.3745
-0.2923
-0.0837
-0.7239
-0.1469
-0.1563
0.2947
-0.0907
-0.1811
0.0301
0.4579
0.3097
0.2775
0.1327
0.2958
0.4443
0.3743
0.1769
0.8769
0.4160
0.3419
-0.0619
0.6534
-0.6500
-0.1013
0.0353
-0.2746
-0.0589
-0.1639
0.7854
-0.4437
-0.4374
0.0429
-0.4342
0.3360
-0.2087
-0.2364
-0.2812
-0.1074
-0.2720
0.3752
0.0039
0.0870
0.0241
-0.0750
0.2155
-0.2159
-0.0187
-0.2019
0.1640
0.1316
-0.4866
-0.1333
0.4158
-0.1819
-0.0130
-0.1416
-0.1798
-0.0139
-0.2570
-0.1240]);


weights = reshape(beta,[8 10]);

% one matrix for each learner (stimulus), for each bootstrap iteration 
% so, not quite as simple as finding a given cell's contribution

figure; 
imagesc(weights')
caxis([-1 1])
cmocean('-curl','pivot',0)
colorbar
title('weight matrix for learner 1, single iteration')
ylabel('Cell #')
xlabel('Stimulus')

print_eps_kp(gcf,fullfile(fn.figs,'ClassAM','AC','Full','SVMweightMatEx_Mar28'));




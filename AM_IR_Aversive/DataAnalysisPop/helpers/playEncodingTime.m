
nlags = 100;

ET = nan(size(GoodData,1),1);
for iUn = 1:size(GoodData,1)
    
halfACF = acf(GoodData(iUn,:)',nlags);
fullACF = [flipud(halfACF); halfACF];
lagvec = [-nlags:-1 1:nlags];


% Fit gaussian 

gaussEqn = 'a*exp(-((x-b)/c)^2)+d';

gaussEqn2 = 'g^2 + p * g * (1/sqrt(4*pi*s^2)) * exp(-(x^2)/(4*s^2))';

f1 = fit(lagvec',fullACF,gaussEqn2,'Normalize','on','Start', [-0.05 -10 50]);
% yeval = f1.g^2 + f1.p * f1.g * (1/sqrt(4*pi*f1.s^2)) * exp(-(lagvec.^2)/(4*f1.s^2));
% 
% figure; 
% % plot(lagvec,fullACF)
% % hold on
% plot(f1,lagvec, fullACF);

% Alt method
% [m,s] = normfit(fullACF);

ET(iUn) = 2*f1.s;

if ET(iUn)>80
    figure;
    plot(f1,lagvec, fullACF);
    title(iUn)
end

end
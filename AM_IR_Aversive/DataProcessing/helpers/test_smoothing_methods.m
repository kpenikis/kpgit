function test_smoothing_methods

InputData = zeros(1,500);
InputData(randperm(500,10))=1;


%~~~~~~~~~~~~~~~~~~
winlen = 500;

figure; hold on
stem(InputData,'k')

for tau = [10 20]
    
    % gaussian first
    clear convwin Output
    
    bs_gaus = tau;
    convwin = gausswin(bs_gaus)';
%     convwin = convwin-min(convwin);
%     convwin = convwin/sum(convwin);
    
    Output = conv(InputData,convwin,'same');
    
    plot(Output,'LineWidth',2)
    
    
    % now exponential
    clear convwin Output
    
    % compute the smoothing constant
    lambda = 2/tau;
    
    %  convwin = exp(-lambda*(linspace(0,1,winlen)));
    convwin = exp(-lambda*(1:winlen));
    %  convwin = convwin./norm(convwin);
    
    Output = conv(InputData,convwin);
    
    plot(Output)
    
%     pct = Output(5+tau)/Output(5)
    
end

 
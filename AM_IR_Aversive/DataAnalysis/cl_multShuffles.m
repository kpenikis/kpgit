
fn = set_paths_directories;
histbinsize = 0.01;


% Reformat shuffled results

Shuffled = struct();

for ii=1:10
    q=load(fullfile(fn.processed,'MPHclassifier',sprintf('ClassData_shuff_%i',ii)));
    DataSh = q.Data; clear q
    
    for iUn=1:size(DataSh,1)
        for irate=1:5
            Shuffled(iUn,irate).dp(:,ii) = DataSh(iUn,irate).Res_L1o.dprime(:,2);
        end
    end
end
clear DataSh


% Judge significance of each MPH pair

q=load(fullfile(fn.processed,'MPHclassifier','ClassData'));
Data = q.Data; clear q

sig_dps = [];

for iUn=1:size(Data,1)
    for irate=1:5
        
        nBeat = sum(Data(iUn,irate).Res_L1o.dprime(:,2) > Shuffled(iUn,irate).dp,2);
        for ib = find(nBeat>9)'
            
            sig_dps = [sig_dps; iUn irate ib Data(iUn,irate).Res_L1o.dprime(ib,2)];
            
        end
        
    end
end

unique(sig_dps(:,1))
unique(sig_dps(:,2))

figure;
plot([0 0],[0 1],'--','Color',0.7.*[1 1 1])
hold on

histogram(sig_dps(:,3),-4:histbinsize:4,'DisplayStyle','stairs',...
    'EdgeColor',0.5.*[1 1 1],'LineWidth',2,'Normalization','cdf');


% Save sig_dps
save(fullfile(fn.processed,'MPHclassifier','sig_dps'),'sig_dps','-v7.3')




% Plot cdf from each uniquely shuffled run of the classifier

figure;
plot([0 0],[0 1],'--','Color',0.7.*[1 1 1])
hold on
for ii=1:10
    
    all_dps = [];
    
    q=load(fullfile(fn.processed,'MPHclassifier',sprintf('ClassData_shuff_%i',ii)));
    DataSh = q.Data; clear q
    
    for iUn=1:size(DataSh,1)
        for irate=1:5
            all_dps = [all_dps; DataSh(iUn,irate).Res_L1o.dprime(:,2)];
        end
    end
    
    ih(ii)=histogram(all_dps,-4:histbinsize:4,'DisplayStyle','stairs',...
    'EdgeColor',0.5.*[1 1 1],'LineWidth',2,'Normalization','cdf');
    
end

xlim([-2 2])
ylim([0 1])




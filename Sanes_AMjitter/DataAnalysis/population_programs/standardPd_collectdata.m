function standardPd_collectdata( rasters, session, cluname, clulabel,...
                        stimpars, behav, indVar )
                    
global InfoTable AmalgamatedSpikes OtherPdSpikes only75

baseFR = calc_baselineFR(rasters);

if only75
    rasters = rasters([rasters.AMdepth]==0.75);
end

these_pds = [2 3 5 6 7];

try
    
for ir = 1:numel(rasters)
    
    % if FR is too low, set data output to nans
    if isempty(rasters(ir).x) || (numel(rasters(ir).x)/max(rasters(ir).y) / (rasters(ir).stimDur/1000)) < 5
        disp('skipping datapoint with too few spikes')
        continue
    end
    
    % Prep for other periodic pd spikes
    rV = jitter_LUT(rasters(ir).AMrate,char(rasters(ir).jitter));
    tV = cumsum( [ ceil(rasters(ir).AMonset+0.75*1000/rV(1)) 1000./rV(2:end) ]);
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    clear x y
    [x,y,prev250FR,prev100FR,prev50FR] = standardPd_getspikes(rasters(ir));
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    for iy = 1:max(y)
        
        this_x = zeros(1,250);
        this_x(x(y==iy)) = 1;
        
        if size(InfoTable,1)+1 > size(AmalgamatedSpikes,1)
            keyboard
        end
        AmalgamatedSpikes(size(InfoTable,1)+1,:) = this_x;
        
        %'Session' 'cluname' 'SU/MU' 'HP' 'LP' 'dBSPL' 'CenterRate' 'BehavState'
        %'depth' 'jitter' 'baselineFR' 'prev250FR' 'prev100FR' 'DiscrmVar' 'trialN'
        DT_addrow = {session cluname ...
            discretize(clulabel, [0.5 1.5 2.5 3.5 4.5], 'categorical',{'unk', 'SU', 'MU', 'noise'})...
            stimpars(1) stimpars(2) stimpars(3) stimpars(4) ...
            categorical(behav, {'P' 'D' 'A'}, {'Passive' 'Drinking' 'Active Behavior'})...
            rasters(ir).AMdepth rasters(ir).jitter baseFR prev250FR prev100FR indVar iy };
        
        if size(InfoTable,2) ~= size(DT_addrow,2)
            keyboard
        end
        
        InfoTable = [InfoTable; DT_addrow];
        
        
        % Also get spikes for other periods of periodic AM
        if strcmp(rasters(ir).jitter,'0')
            for ipd = these_pds
                clear x_tr x_idx
                x_tr = rasters(ir).x(rasters(ir).y==iy);
                x_idx = x_tr(x_tr>=tV(ipd-1) & x_tr<=tV(ipd)) - tV(ipd-1) + 1;
                OtherPdSpikes( size(InfoTable,1), x_idx , these_pds==ipd ) = 1;
            end
        end
        
    end
end

catch
    keyboard
end


end




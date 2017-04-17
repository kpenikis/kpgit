function N = get_drinking_trials(block)
    

%  [ block n_drinking_trs ] 
N_drinking_trials = [...
        17 21;
        21 24;
        22 97;
        23 0;
        27 106;
        31 117;
        35 109;
        36 0;
        40 107;
        41 0;
        45 124;
        46 0;
        50 139;
        55 216;
        61 0;
        62 0;
        63 0;
        64 148;
        68 118;
        69 0;
        73 193;
        79 31;
        80 0;
        88 153;
        89 0;
        90 0;
        98 150;
        99 0;
       100 0;
       106 75;
       107 0;
       108 0;
       112 0;
       116 0;
       122 41;
       128 56;
       133 144;
       138 nan;
       139 86;
       141 198;
       146 0;
       147 0;
       148 57;
       157 101;
       158 52;
       159 0;
       165 197;
       169 214;
       174 147;
       178 159  ];

if sum(N_drinking_trials(:,1)==block)<1
    disp('block not found in vector')
    keyboard
end

N = N_drinking_trials(N_drinking_trials(:,1)==block,2);

end
 



function checkSortingQuality( BLOCK, subject, session_label, tetN )
%    IN PROGRESS...
%  KP, 2017-12
%


global fn

fn = set_paths_directories(subject,session_label);


%% Get raw physiology data from epData and filter

Phys = get_phys_data(BLOCK,subject,tetN);



%% Plot Wave (also plot pre filtered to create way to remove noise in
% beginning of session)

% 




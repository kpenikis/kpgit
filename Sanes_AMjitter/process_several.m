function process_several

blocks      = {  64   [68:69]  [141:142] [178:181] };
sess_labels = { 'JA'   'JC'      'MG'      'PA'    };

for isess = 1:numel(sess_labels)
    
    try
        close all
        pp_prepare_format(blocks{isess},'IIIf_230115',sess_labels{isess})
        
    catch
        keyboard
    end
    
end


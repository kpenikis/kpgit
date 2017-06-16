function process_several

blocks      = { 157   [122:125]  [148]   [165:168]  [174:177] };
sess_labels = { 'OB'   'MC'       'NA'      'OE'       'OH'   };

for isess = 1:numel(sess_labels)
    
    try
        close all
        pp_prepare_format(blocks{isess},'IIIf_230115',sess_labels{isess})
        
    catch
        keyboard
    end
    
end


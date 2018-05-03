function print_svg_kp(fig_handle,savelocation,YESpleasePRINT)

if nargin < 3
    YESpleasePRINT = 1;
end

if YESpleasePRINT
    set(figure(fig_handle),'PaperPositionMode','auto');
    print(figure(fig_handle),'-painters','-dsvg', savelocation)
end

end
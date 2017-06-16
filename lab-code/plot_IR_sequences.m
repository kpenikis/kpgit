function plot_IR_sequences

load IRsequences.mat

colors = hsv(4);

xvec = [ 0 1 1 2 2 3 3 4 4 5 5 6 ];

hf(1)=figure;
plot(xvec,reshape([IR_A;IR_A],1,12),'Color',colors(1,:),'LineWidth',10)
set(gca,'yscale','log','xtick',[],'ytick',[])
xlim([0 6])
print(hf(1),'IR_A_plotted','-dpdf','-bestfit')

hf(2)=figure;
plot(xvec,reshape([IR_B;IR_B],1,12),'Color',colors(2,:),'LineWidth',10)
set(gca,'yscale','log','xtick',[],'ytick',[])
xlim([0 6])
print(hf(2),'IR_B_plotted','-dpdf','-bestfit')

hf(3)=figure;
plot(xvec,reshape([IR_C;IR_C],1,12),'Color',colors(3,:),'LineWidth',10)
set(gca,'yscale','log','xtick',[],'ytick',[])
xlim([0 6])
print(hf(3),'IR_C_plotted','-dpdf','-bestfit')

hf(4)=figure;
plot(xvec,reshape([IR_D;IR_D],1,12),'Color',colors(4,:),'LineWidth',10)
set(gca,'yscale','log','xtick',[],'ytick',[])
print(hf(4),'IR_D_plotted','-dpdf','-bestfit')



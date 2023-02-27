close all;
clear all;
clc;

load radiation.mat

E = sqrt(abs(E4_x).^2 + abs(E4_y).^2 + abs(E4_z).^2);
I = E.^2;

figure(2)
image(yy*2.38E-3,xx*2.38E-3,I./max(max(I)),'CDataMapping','scaled')
colorbar;
axis square;
sum(E.^2,'all') 
% figure(2)
% image(z,x,real(EE),'CDataMapping','scaled')
% axis square;

xlabel('x (mm)','fontsize',50);
ylabel('y (mm)','fontsize',50);
% xlim([5,25]); % adjust as needed
% ylim([-15,20]); % adjust as needed
% title(['\theta = ',num2str(0.71),' deg']);
% set(gca,'xtick',0:2:20);
% set(gca,'ytick',0:5:30);
set(gca,'fontsize',50,'fontname','Arial','linewidth',2)
% set(gca,'position',[0.1,0.1,0.8,0.4])
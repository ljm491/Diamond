clear all;
close all;
clc;

PM_number = 1;

Color(1,:) = [46 47 138]./255;
Color(2,:) = [118 105 175]./255;
Color(3,:) = [196 174 208]./255;
Color(4,:) = [210 108 129]./255;
Color(5,:) = [199 32 38]./255;

ClearColor(1,:) = [109 110 173]./255;
ClearColor(2,:) = [173 165 207]./255;
ClearColor(3,:) = [240 226 247]./255;
ClearColor(4,:) = [233 182 192]./255;
ClearColor(5,:) = [210 120 120]./255;

c = 2.99792458;

diamond_tilt = 45 * pi / 180;

q(:) = [0.418 0.5095 0.585 0.6495 0.706];
q1 =q(PM_number) * pi / 180; % rad

f4 = 2.5:0.1:25; % THz

d = 500 / 1000 / 10; % cm

w4 = 104; % um

k1=n_diamond(206) * 2 * pi * 206E+12 /c/10^10; %cm-1
k2=n_diamond(166) * 2 * pi * 166E+12 /c/10^10;

k3 =n_diamond(f4+40) .* 2 .* pi .* (f4+40).*1E+12 ./c./10^10;
k4 =(n_diamond(f4)-0.0) .* 2 .* pi .* f4.*1E+12 ./c./10^10;

idl_angle = asin(sin(diamond_tilt)/n_diamond(166));

q2 = pi - idl_angle - q1;
l2 = sin(idl_angle) * (k2 + k3) / sin(q2);
q3 = asin( sin(idl_angle + q1) .* (k1-l2) ./ k4);

qd = idl_angle + q3;
qd = qd - idl_angle;
qe = asin (sin(qd).*n_diamond(f4));

qsp = (c.*100./f4) ./ pi ./ n_diamond(f4) ./ w4;
qsp_air =  (c.*100./f4) ./ pi ./ w4;

q_shade_x = [f4(1:length(f4)) f4(length(f4):-1:1)];
q_shade_y = [qd+qsp qd(length(qd):-1:1)-qsp(length(qd):-1:1)];

% TIR_crit_angle = asin(1./n_diamond(f4));
% Brewster_angle = atan(1./n_diamond(f4));

figure(1)

fill(q_shade_x,q_shade_y.*180./pi,'r','FaceColor',ClearColor(PM_number,:),'EdgeColor','None');hold on;
plot(f4,qd.*180./pi,'-','linewidth',4,'Color',Color(PM_number,:));hold on;
% plot(f4,TIR_crit_angle.*180./pi,'--','color',[0.4 0.4 0.4],'linewidth',4);
% plot(f4,Brewster_angle.*180./pi,'--','color',[0.7 0.7 0.7],'linewidth',4);

xlabel('Frequency (THz)','fontsize',50);
ylabel('Internal angle (deg)','fontsize',50);
xlim([5,25]); % adjust as needed
ylim([-15,20]); % adjust as needed
title(['\theta = ',num2str(0.71),' deg']);
% set(gca,'xtick',0:2:20);
% set(gca,'ytick',0:5:30);
set(gca,'fontsize',50,'fontname','Arial','linewidth',2)
set(gca,'position',[0.1,0.3,0.8,0.4])
set(gca,'linewidth',4)
% grid on

% DR = load('Ang Spread 0.418 THz.txt');
% f4_DR = DR(:,1);
% qd_DR = DR(:,2)-idl_angle/pi*180;
% qsp_DR = DR(:,3);
% 
% scatter(f4_DR,qd_DR,200,Color(PM_number,:),'filled');
% Err = errorbar(f4_DR,qd_DR,qsp_DR,qsp_DR,'o','MarkerSize',0.1);
% Err.Color = Color(PM_number,:);
% Err.LineWidth = 4;
% Err.CapSize = 20;
% 
% legend('Spread (PWA)','Central (PWA)','Central (NS)','Spread (NS)','location','best');

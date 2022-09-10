clear all;
close all;
clc;

%% Set parameters
PM_number_list = 1:1:5;

Curve_type = '-';
Color(1,:) = [46 47 138]./255;
Color(2,:) = [118 105 175]./255;
Color(3,:) = [196 174 208]./255;
Color(4,:) = [210 108 129]./255;
Color(5,:) = [199 32 38]./255;

for PM_number = PM_number_list

d = 500 / 1000 / 10; % diamond thickness/cm

q(:) = [0.418 0.5095 0.585 0.6495 0.706];

q1 =q(PM_number) * pi / 180;

diamond_tilt = 45 * pi / 180;

f_MIR(:) = [5 7.5 10 12.5 16];

MIR_center = 1332.5 + 33.3 * f_MIR(PM_number);
MIR_bdw = 33.3 * 6;
MIR_energy = 8.4;

%%

f4 = 2.5:0.1:30;

c = 2.99792458;

k1 = n_diamond(206)   * 2  * pi  * 206E+12         /c /10^10;
k2 = n_diamond(166)    * 2  * pi  * 166E+12         /c /10^10;
k3 = n_diamond(f4+40) .* 2 .* pi .* (f4+40).*1E+12 ./c./10^10;
k4 = n_diamond(f4)    .* 2 .* pi .* f4     .*1E+12 ./c./10^10;

idl_angle = asin(sin(diamond_tilt)/n_diamond(166));

q2 = pi - idl_angle - q1;
l2 = sin(idl_angle) * (k2 + k3) / sin(q2);
q3 = asin( sin(idl_angle + q1) .* (k1-l2) ./ k4);

qd = idl_angle + q3;

k23_z = (k2 + k3) .* sin(q1) ./ sin(pi-idl_angle-q1);
k14_z = k4 .* sin(idl_angle + q1 + q3) ./ sin(q2);

dk = k23_z - k14_z; % dk_x = 0;

pm = sin (dk.*d./2) .^2 ./ (dk.*d./2).^2; % sinc function squared

figure(1)
plot(f4,pm,Curve_type,'linewidth',4,'color',Color(PM_number,:)); hold on;
set(gca,'fontsize',30,'fontname','Arial','linewidth',2)
% title('0.5-mm-thick diamond','fontsize',30);
xlabel('Frequency (THz)','fontsize',30);
ylabel('sinc^2(\Deltak_zL/2)','fontsize',30);
xlim([2.5,30]); % adjust as needed
ylim([0,1]); % adjust as needed
legend({'0.42','0.51','0.59','0.65','0.71'},'fontsize',30,'Location','best')
set(gca,'position',[0.1,0.1,0.8,0.7])

qd = qd - idl_angle;

qe = asin (sin(qd).*n_diamond(f4));

figure(2)
plot(f4,qe.*180./pi+diamond_tilt*180/pi-diamond_tilt*180/pi,Curve_type,'linewidth',6,'color',Color(PM_number,:));hold on;
xlabel('Frequency (THz)','fontsize',50);
ylabel('External angle (deg)','fontsize',50);
xlim([2,20]); % adjust as needed
ylim([-32,48]); 
% set(gca,'xtick',0:2:20);
% set(gca,'ytick',0:15:100);
set(gca,'fontsize',30,'fontname','Arial','linewidth',2)
legend({'0.42','0.51','0.59','0.65','0.71'},'fontsize',30,'Location','best')
set(gca,'position',[0.1,0.1,0.8,0.35])

I_MIR_shifted = MIR_energy * exp( -2.*((f4 -(MIR_center-1332.5)/33.33) ./(MIR_bdw/33.33)).^2) /MIR_bdw;
I_MIR_shifted = I_MIR_shifted./max(I_MIR_shifted);
I_THz = I_MIR_shifted .* pm .* f4.^2;

I_THz = I_THz .* cos(qd).^2; % projection of P

reflec = (n_diamond(f4) .* cos(qe) - cos(qd)) ./ (cos(qd) + n_diamond(f4) .* cos(qe));
I_THz = I_THz .* (1-reflec.*conj(reflec)); % Fresnel reflection loss

I_THz = I_THz ./ max (I_THz);

% I_THz = I_THz ./ 0.8904 *2.8/10000 /1.06^6;

% chi3 = sqrt(0.0015/real(sum(I_THz,'all')/sum(I_MIR_shifted,'all')));

% dlmwrite('output.txt',[f4' I_THz'],'\t');

% ref = load('10THz_DRsimu.txt');
% ref_freq = ref(:,1);
% ref_I = ref(:,2);
% ref_I = ref_I./max(ref_I);
% ref_I=ref_I.^2.*0.35;

figure(4)
plot(f4,I_THz,Curve_type,'linewidth',4,'color',Color(PM_number,:));hold on;
% plot(ref_freq,ref_I,'bo','linewidth',2);
set(gca,'fontsize',30,'fontname','Arial','linewidth',2)
xlabel('Frequency (THz)','fontsize',30);
ylabel('Intensity (a.u.)','fontsize',30);
xlim([min(f4),25]); 
% ylim([0,1]);
% grid on;
legend({'0.42','0.51','0.59','0.65','0.71'},'fontsize',30,'Location','best')
set(gca,'position',[0.1,0.1,0.8,0.8])

end

% 
% Calculation of nonlinear polarization P(r,omega) of FWM in diamond
% 
% Based on Eq. 8.6 of Principles of Nonlinear Optics by Y.R.Shen
%
%
% Input beam symbolized by 
% 1: Signal     1.45 um     narrowband
% 2: Idler      1.8 um      narrowband
% 3: MIR        6 um        broadband & tunable 
% 4: THz
% 
% Diamond coordinate:
% x: parrrel to TOPAS output at exit; points opposite to TOPAS output
% y: vertical; points up
% z: propagation direction of idler and MIR
%
% P(3)_x = chi(3)_xxxx E1_x E2_x E3_x
% 

close all;
clear all;
clc;

c = 2.99792458;

Diamond_dimension_transverse = 500; % um; > input beam diameter
Diamond_thickness = 500; % um

d_r_transverse = 10; % um; << input beam diameter
d_z = 5; % um; << THz wavelength;

x = (-Diamond_dimension_transverse/2:d_r_transverse:Diamond_dimension_transverse/2)';
y = -Diamond_dimension_transverse/2:d_r_transverse:Diamond_dimension_transverse/2;
z = -Diamond_thickness/2:d_z:Diamond_thickness/2;

E1_0 = 1; % amplitude
E2_0 = 1;
E3_0 = 1;

Theta = 0.585 * pi / 180;

f1 = 206; % THz
k1 = (n_diamond(f1))*2*pi/(c * 100 / f1); % per um
w1 = 170; % um
zR1 = pi * w1 ^ 2 * (n_diamond(f1)) / (c * 100 / f1); % um

f2 = 166; % THz
k2 = n_diamond(f2)*2*pi/(c * 100 / f2); % per um
w2 = 210; % um
zR2 = pi * w2 ^ 2 * n_diamond(f2) / (c * 100 / f2); % um

f3 = 50; % THz
k3 = n_diamond(f3)*2*pi/(c * 100 / f3); % per um
w3 = 170; % um
zR3 = pi * w3 ^ 2 * n_diamond(f3) / (c * 100 / f3); % um

for jj = (1:length(x))
    for ii = (1:length(z))
        E1(jj,:,ii) = cos(Theta) * E1_0 /(1 + i * z(ii) * cos(Theta) / zR1 - i * x(jj) * sin(Theta) / zR1) * exp(-((x(jj)*cos(Theta)+z(ii)*sin(Theta))^2+y.^2)./w1^2./(1 + i * z(ii) * cos(Theta) / zR1 - i * x(jj) * sin(Theta) / zR1)) * exp(i*k1*z(ii)*cos(Theta) - i*k1*x(jj)*sin(Theta));
    end
end

for ii = (1:length(z))    
    E2(:,:,ii) = E2_0 / (1 + i * z(ii) / zR2) * exp(-((x).^2 + y.^2)./w2^2 ./(1 + i * z(ii) / zR2)) * exp(i*k2*z(ii));
    E3(:,:,ii) = E3_0 / (1 + i * z(ii) / zR3) * exp(-((x).^2 + y.^2)./w3^2 ./(1 + i * z(ii) / zR3)) * exp(i*k3*z(ii));
end

P =conj(E1).*E2.*E3;

% clear E1 E2 E3


save polarization
Plot_Polarization

%close all;

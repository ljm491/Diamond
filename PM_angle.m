clear all;
close all;
clc;

c = 2.99792458;

f4 = 2.5:0.1:30; % set phase matched THz frequency

k1=n_diamond(206) * 2 * pi * 206E+12 /c /10^10;
k2=n_diamond(166) * 2 * pi * 166E+12 /c /10^10;

k3 =n_diamond(f4+40) * 2 * pi .* (f4+40) * 1E+12 / c /10^10;
k4 =n_diamond(f4)    * 2 * pi .* f4      * 1E+12 / c /10^10;

q32 = 0 * pi / 180;

k0 = sqrt(k2^2 + k3.^2 - 2* k2 * k3 * cos(pi-q32) );

q02 = acos( (k2^2 + k0.^2 - k3.^2)/2/k2./k0 );

q01 = acos( (k0.^2 + k1^2 - k4.^2)/2/k1./k0 );

q12 = q01 + q02;

q40 = acos((k0.^2+k4.^2-k1^2)/2./k4./k0);

fprintf('Theta = %f\n\n',(q01 + q02) * 180 / pi)
fprintf('Alpha = %f\n\n',(q40 + q02) *180 / pi)

plot(f4,(q01 + q02) * 180 / pi)

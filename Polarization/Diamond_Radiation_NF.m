close all;
clear all;
clc;

load polarization

Det_Plane_x_range = 1 * 1000; % um

Det_Plane_y_range = 1 * 1000; % um
Det_Plane_y_center= 0*1000;

Det_Plane_dr = 1/100 * 1000; % um

xx = -Det_Plane_x_range/2 : Det_Plane_dr : Det_Plane_x_range/2;
yy = Det_Plane_y_center-Det_Plane_y_range/2 : Det_Plane_dr :Det_Plane_y_center+Det_Plane_y_range/2;
zz = 0 + 0.750 * 1000; % um

[Xr2,Yr2,Zr2] = meshgrid(x,y,z);
f4 = f3 - (f1 - f2); % THz
k4 = n_diamond(f4)*2*pi/(c * 100 / f4); % per um

E4_x = zeros(length(xx),length(yy));
E4_y = zeros(length(xx),length(yy));
E4_z = zeros(length(xx),length(yy));
for iii = (1:length(xx))
    for jjj = (1:length(yy))

        xxi = xx(iii);
        yyj = yy(jjj);
        r = sqrt(xxi^2+yyj^2+zz^2);
        dismat = sqrt((xxi-Xr2).^2+(yyj-Yr2).^2+(Zr2-zz).^2);
        Phase_factor = exp(1i*k4*dismat)./dismat;
        E4_x(iii,jjj) = sum(P.*(1-xxi^2/r^2).*Phase_factor,'all');
        E4_y(iii,jjj) = sum(P.*(-xxi*yyj/r^2).*Phase_factor,'all');
        E4_z(iii,jjj) = sum(P.*(-xxi*zz/r^2).*Phase_factor,'all');
    end
end

save radiation
Plot_Radiation


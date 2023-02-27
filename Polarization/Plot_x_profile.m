close all;
% clear all;
clc;

load radiation

Ix = E(:,26).^2;
Iy = E(26,:).^2
Ix = Ix./max(Ix);
Iy = Iy./max(Iy);
plot(yy,Ix)
hold on
plot(xx,Iy)
close all;
clear all;
clc;

load polarization

my = round(length(y)/2);

for jj = (1:length(x))
    for ii = (1:length(z))
        EE(jj,ii)=P(jj,my,ii);
    end
end

figure(1)
image(z,x,abs(EE),'CDataMapping','scaled')
axis square;
figure(2)
image(z,x,real(EE),'CDataMapping','scaled')
axis square;


for jj = (1:length(x))
    for ii = (1:length(y))
        EEE(jj,ii)=P(jj,ii,100);
    end
end
figure(3)
image(y,x,abs(EEE),'CDataMapping','scaled')
axis square;
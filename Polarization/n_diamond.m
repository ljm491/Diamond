function [n] = n_diamond(nu)
% return ref. index
% nu frequency in the unit of THz
l = 299.792458./nu.*1000; % unit: micron
n = sqrt(1+(0.3306*l.^2)./(l.^2-175^2)+(4.3356.*l.^2)./(l.^2-106^2));
%n = 2.3804;
end


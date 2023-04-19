function [ht, t] = r_rcro(kt, Tb, spb, r, ts)
 if ~exist('ts','var')
     % third parameter does not exist, so default it to something
      ts = 0;
 end
t = -kt*Tb+ts:Tb/spb:kt*Tb+ts-Tb/spb;
ht = zeros(size(t));
R = 1/Tb;

% get rid of floating point errors
t_tb_4r = abs(abs(t)-Tb/(4*r)) < 1e6*eps(min(abs(t), Tb/(4*r)));

% PW func
ht(t==0) = 1-r+4*r/pi;
ht(t_tb_4r) = r/sqrt(2)*((1+2/pi)*sin(pi/(4*r)) + (1-2/pi)*cos(pi/(4*r)));
t(abs(t)==Tb/(4*r));

cond = t~=0 & ~t_tb_4r;
ht(cond) = (sin(pi.*R.*t(cond).*(1-r)) + 4.*R.*r.*t(cond).*cos(pi.*R.*t(cond).*(1+r)))./(pi.*R.*t(cond).*(1 - (4.*R.*r.*t(cond)).^2));


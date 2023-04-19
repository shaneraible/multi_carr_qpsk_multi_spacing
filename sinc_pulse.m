function [p, t] = sinc_pulse(Tb, n, spd)

    t =  -n*Tb:Tb/spd:n*Tb;
    R = 1/Tb;
    p = sin(pi*R*t)./(pi*R*t);
    p(t==0)=1;
end
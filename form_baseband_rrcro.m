function [signal, t, ak, bits] = form_baseband_rrcro(N, R, spb, kt, r, A)
    ak = rand_bits(N);
    ak(ak==0) = -1;
    ak = ak*A;
    start = 1;
    
    Tb  = 1/R;

    signal = zeros(1,N*spb+kt*spb*2 - spb);
    bits = zeros(1,N*spb+kt*spb*2 - spb);
    t = -kt*Tb:Tb/spb:(Tb*N+kt*Tb*2 - Tb/spb - kt*Tb - Tb);
    % figure()
    
    for k = 1:1:N
        last = kt*spb*2+start-1;
        
        curr_pulse = ak(k) * r_rcro(kt, Tb, spb, r);
        bits(bits==0)=nan;
        bits(start+(last-start+1)/2 - 1) = ak(k);
        
        curr_signal = signal(start:last);
        signal(start:last) = curr_pulse+signal(start:last);
        start = start + spb;
    end

end
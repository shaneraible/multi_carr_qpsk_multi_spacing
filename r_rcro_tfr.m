function hf = r_rcro_tfr(f, B, r)


f0 = B/(r+1);
fd = f0 * r;
f1 = f0-fd;

hf = zeros(1, length(f));
eq_tmp = (abs(f)<B) & (abs(f)>=f1);

hf(abs(f)<f1)               = 1;
hf(eq_tmp)  = .5 * (1 + cos((pi*(abs(f(eq_tmp)) - f1))/(2*fd)));
hf(abs(f)>=B)               = 0;

end

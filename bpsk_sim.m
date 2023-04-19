
close all
clear

%% B %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N   = 89;       % num bits 
R   = 1e6;     % bit rate    
spb = 16;       % symbols per bit
kt  = 20;        % tail size
r   = 0.5;      % rolloff factor
A_bb = 1;       % baseband amplitude
Tb  = 1/R;      % bit time
max_noise_var = 12;


fs = (1/(Tb/spb));
fc = fs/4;       % carrier frequency
wc = 2*pi*fc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[m_t, t, ak, bits] = form_baseband_rrcro(N, R, spb, kt, r, A_bb);

figure()
subplot(2,1,1)
stem(1:N, ak)
title("Random bits, ak")
xlabel("k")
ylabel("ak")
subplot(2,1,2)
plot(t, m_t)
hold on
plot(t, bits, 'x');
title("Baseband Signal, R-RCRO, r=.5")
xlabel("t");
ylabel("m(t)");
hold off

i_t  = Tb/spb; % multiplier for time to index
t_shift_idx = spb;

rrcro_conv = r_rcro(kt, Tb, spb, r, 0);
rrcro_conv = rrcro_conv(end:-1:1);

k   = 1;
n0  =2;
c   = 2*k/n0;
yt  = conv(m_t, rrcro_conv, 'same');

shift = 0;
figure()
plot(t, yt)
hold on
plot(t, bits.*abs(yt), 'x');
title('Match-filtered R-RCRO Baseband Signal')
xlabel("t");
ylabel("m(t)");

figure()
subplot(2,3,1)
stem(t, bits)
xlabel('t')
ylabel('at')
title('Bits transmitted')
noise = randn(size(m_t));

subplot(2,3,2)
hold off
plot(t, m_t)
xlabel('t')
ylabel('m(t)')
title('Baseband Signal')
    
cos_func = cos(wc.*t + pi/2);

for i=0:1:max_noise_var
    noise_amplitude = i;

    s_t_bpsk = m_t.*cos_func;
    
    subplot(2,3,3)
    hold off
    plot(t, s_t_bpsk)
    xlabel('t')
    ylabel('s(t)')
    title('BPSK Modulated Signal')
    
    r_t = s_t_bpsk + noise_amplitude*noise;
    
    subplot(2,3,4)
    hold off
    plot(t, r_t)
    xlabel('t')
    ylabel('r(t)')
    title(['Received Signal -- noise σ = ' num2str(noise_amplitude)])
    
    r_t_baseband = r_t .* cos_func;
    r_t_baseband = lowpass(r_t_baseband,2*R,fs);

    subplot(2,3,5)
    plot(t, r_t_baseband)
    xlabel('t')
    ylabel('r(t)')
    title(['Baseband Rx Signal -- σ = ' num2str(noise_amplitude)])

    subplot(2,3,6)
    hold off
    yt_noise = conv(r_t_baseband, rrcro_conv, 'same');
    bits_recovered = (abs(bits).*yt_noise>0) + -1*(abs(bits).*yt_noise<0);
    bits_recovered(abs(bits_recovered)~=1) = nan;
    bit_errors = sum(bits ~= bits_recovered & ~isnan(bits_recovered));
    
    % locations of bit errors
    error_bits = double(((bits_recovered~=bits) & ~isnan(bits_recovered)));
    error_locs = (abs(error_bits)==0);
    error_bits(abs(error_bits)==0) = nan;
    
    plot(t,yt_noise)
    xlabel('t')
    ylabel('y(t)')
    title(['Filtered Received Signal -- with noise σ = ' num2str(noise_amplitude)])
    hold on
    plot(t, (bits_recovered).*abs(yt_noise), 'g*');
    hold on
    plot(t, error_bits.*(yt_noise), 'r*');

    if noise_amplitude + 1 <= max_noise_var
        disp(['Press any key to proceed to σ = ' num2str(noise_amplitude+1) ' of ' num2str(max_noise_var)]);
    else
        disp('Press any key to continue...')
    end
    pause;
    hold off
end

%% plotting the PSD and Spectrum
% N = 89;
% 
% f0  = 1e6;  % bit rate
% spb = 16;   % samples per bit
% kt  = 20;   % bit times
% r   = .5;   % rolloff factor
% ts  = 0;    % time shift
% Tb  = 1/f0; % bit time
% A   = 10;

% (40*1/1e6+89*1/1e6-1/1e6)/(1/1e6/16)
f0 = R;
fs = f0*spb;
% fc = fs/4;       % carrier frequency
wc = 2*pi*fc;
B = f0 * (1+r)/2;

psf_avg = zeros(1, (2*kt*Tb+N*Tb-Tb)/(Tb/spb));
figure()
disp('running PSD simulation...')

for p = 0:3
    num_trials = 10^(p+1);

    for i = 1:num_trials
        [m_t, t, ak, bits] = form_baseband_rrcro(N, f0, spb, kt, r, A_bb);
        cos_func = cos(wc.*t + pi/2);

        s_t = m_t.*cos_func;
        Sf = fft(s_t);
        Sf = Sf * Tb/spb;
        Sf = fftshift(Sf);    
        ls = length(Sf);
        fs = f0*spb;
        f = (-ls/2:ls/2-1)/ls*fs;
        
        psf = abs(Sf).^2 ./ (N * Tb);
        psf_avg = psf_avg + psf;
    end
    psf_avg = psf_avg / num_trials;
    subplot(2,2, p+1)
    plot(f, psf_avg)
    title(['BPSK w/ R-RCRO Pulse PSD, ' num2str(num_trials) ' trials'])
    xlabel('f (Hz)')
    ylabel('Ps(f)')
    
    psf_theor = .25 * A_bb^2 * Tb * (r_rcro_tfr(f-fc, B, r) + r_rcro_tfr(f+fc, B, r));
    hold on
    plot(f, psf_theor)
    hold off
    legend('Simulation','Theoretical')

end




%% plotting bit error rates 
% N = 1000;

start_ebno = -40;
end_ebno   = 20;
delta_ebno = .1;
num_trials  = (end_ebno-start_ebno)/delta_ebno + 1;

ebn0_trials = zeros(1, num_trials);
ber_trials  = zeros(1, num_trials);

avg_trials = 50;
fs = (1/(Tb/spb));
eb = Tb*A_bb^2 *  (1/2);
ebno_list = start_ebno:delta_ebno:end_ebno;
sigma_list = sqrt(eb./(10.^(ebno_list./10)) .* fs ./ 2);


disp('running eb/n0 simulation...')
tic
for t_idx = 1:avg_trials
    noise = randn(size(m_t));
    trial_num = 1;
    [m_t, t, ak, bits] = form_baseband_rrcro(N, R, spb, kt, r, A_bb);
    cos_func = cos(wc.*t + pi/2);
    samps = abs(bits).*conv(m_t, rrcro_conv, 'same');
    s_t = m_t.*cos_func;

    for ebno = ebno_list

        sigma = sigma_list(trial_num);

        rx_signal = s_t + sigma*noise;
        rx_basband = rx_signal.*cos_func;

        yt_noise = conv(rx_basband, rrcro_conv, 'same');
        
        bits_recovered = (abs(bits).*yt_noise>0) + -1*(abs(bits).*yt_noise<0);
        bits_recovered(abs(bits_recovered)~=1) = nan;
        bit_errors = sum(bits ~= bits_recovered & ~isnan(bits_recovered));
        
        ber = bit_errors / N;
    
        ebn0_trials(trial_num) = ebno;
        ber_trials(trial_num) = ber_trials(trial_num) + ber;
        trial_num = trial_num + 1;
    end

end
time = toc;

figure()

ber_trials = ber_trials/avg_trials;

% expected_ber = qfunc(sqrt(2*10.^(ebn0_trials/10)));
expected_ber = qfunc(sqrt(2*10.^(ebn0_trials/10)));
hot
plot(ebn0_trials, ber_trials)
hold on
plot(ebn0_trials, expected_ber)
legend('Observed','Expected')
xlabel('Eb/N0')
ylabel('Pe')
title(['Pe vs Eb/No for ' int2str(avg_trials) ' trials with ' int2str(N) ' points in ' int2str(time) 's'])



% now we create a header
Hd= [1 -1 1 -1 -1 1 1 1 -1 -1 -1 1 -1 1 1 -1 -1 -1 1 1];% header bits

%Hd= [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];% header bits
i_p_order = 3;
fs_h=60e6; % sampling frequency
ts_h=1/fs_h; % sampling time 
fc_h=20e6; % bandwidth (same as of the signal)
tc_h=1/fc_h; % bandwdith rate
samples_symbol_header=fs_h/fc_h;
t_h=0:ts_h:tc_h-ts_h;
h_pulse=cos(2*pi*fc_h*t_h);

t_h = (0:1:(length(Hd)*samples_symbol_header-1))*ts_h; % Total time by which N_symbol send

header_modulate_signal = zeros(1,length(t_h)); % Storage for generation of Modulated signal

for ik=1:length(Hd)
    header_modulate_signal((ik-1)*samples_symbol_header+1:ik*samples_symbol_header)=Hd(ik)*h_pulse;
end

% to use with the OFDM_signal: 
OFDM_Interpolated = interp(OFDM,i_p_order);

Transmitted_signal = [header_modulate_signal OFDM_Interpolated];
t_new_2 = (0:1:(length(Transmitted_signal)-1))*ts;
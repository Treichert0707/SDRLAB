
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% removing the header in the receiver

%assuming OFDM_received_signal is the received signal after the Pluto SDR 
amplitude = OFDM_received_signal;
y=xcorr(header_modulate_signal,amplitude); % correlation function
[m,ind]=max(y);

index_time_start=length(amplitude)-ind+(length(Hd)*samples_symbol_header);
index_time_end = index_time_start+(3*(size(x_modulating_with_cpe,1)*size(x_modulating_with_cpe,2)));

%amp_info=amplitude(index_time_start+1:index_time_end);
amp_info=amplitude(index_time_start+1:index_time_end);
N1=length(amp_info);
t = (0:1:(length(amp_info)-1))*ts; % Total time by which N_symbol send
xn_channel=amp_info;
%% removing the interpolated data 
OFDM_received_2 = xn_channel(1);
 for jj=1:((length(xn_channel)-3)/3) %interpolation order
     OFDM_received_2=[OFDM_received_2 xn_channel(3*jj)];
 end


%% Parameters
load('molecular_absorption.mat')
tx_gain = 26; % TX gain in dB
rx_gain = 26; %RX gain in dB
k_b = 1.38064852*10.^(-23); %Boltzmann's Constant 
T = 280; % temperature in Kelvin 
B = 10.^(9); % Bandwidth of the signal 
noise_power = 1.9*10.^(-17)*B; % Power of the noise (psd = 1.9.*10^(-17) W/Hz)
power_tx = 100.*10.^(-3); %Transmitted power in W 
conversion_loss = 7; % 7 dB conversion loss 
noise_factor = 8.5; %8.5 dB NF 
%%
r = [0.1:0.1:2]; % link distance from 10cm to 2 metres 

freq_2_use = f(1312); %f('value to match with 1 THz') (Look within the values from the absorption.mat file)
k_2_use = k(1312); %matching k value. 
% find the path loss 
for i = 1:length(r)
    for j = 1:length(freq_2_use)
    ab_loss(i,j) = 1./(exp(-k_2_use(j).*r(i)));
    ab_loss_db (i,j) = 10*log10(ab_loss(i,j));
    end
end
c_const = 3e8; % speed of light 
wavelength = c_const./freq_2_use ; % find the wavelength from the frequencies.. 
A_eff = wavelength.^2./(4.*pi); % effective aperture of a dipole
 for i = 1:length(r)
     for j = 1:length(A_eff)
         spr_loss(i,j) = 4.*pi.*r(i).^2./(A_eff(j));
         spr_loss_db(i,j) = 10*log10(spr_loss(i,j)); % find the spreading loss in dB
         path_loss_db(i,j) = spr_loss_db(i,j) + ab_loss_db(i,j); % total path loss in dB
     end
 end
Power_tx_dbw = 10*log10(power_tx);
for i = 1:length(path_loss_db)
   total_loss(i) = tx_gain + rx_gain - path_loss_db(i);
   power_rx(i) =10^((Power_tx_dbw+total_loss(i))/10); %received power (Linear value)
   power_rx_dbm(i) = 10*log10(power_rx(i)/0.001); % received power in dB and convert to dBm
   snr(i) = power_rx(i)./noise_power; % linear SNR 
   snr_db(i) =10*log10(snr(i)); % SNR in dB
   power_if(i) = power_rx(i)./(10.^(conversion_loss/10)); % convert 7 dB to linear and divide in the received power
   power_if_db(i) = 10.*log10(power_if(i)./(0.001)); % make it dB and convert to dbM
   snr_if_db(i) = snr_db(i)-noise_factor; % in dB, subtract the noise factor from the receiver SNR 
   snr_if(i) = 10.^(snr_if_db(i)/10); %convert the dB value to a linear value 
end


%%
figure(6)
plot(r,power_rx_dbm, 'Linewidth', 3.5)
title('Received power')
xlabel('distance (m)')
ylabel('Receiver power (dBm)')
%%
figure(7)
plot(r,power_if_db, 'Linewidth', 3.5);
title('Received power after the mixer')
xlabel('distance (m)')
ylabel('Receiver power (dB)')
%%
figure(8)
plot(r,snr_db, 'Linewidth', 3.5)
title('SNR')
xlabel('distance (m)')
ylabel('SNR at the receiver (dB)')

%%
figure(9)
plot(r, snr_if_db, 'Linewidth', 3.5)
title('SNR')
xlabel('distance (m)')
ylabel('SNR at the receiver after the mixer(dB)')

%%
    data_bpsk = B; % data rate in bits per second

for i= 1:length(snr)
   Eb_by_N0_bpsk(i) = snr_if(i)*B/(data_bpsk); % the linear expression (use the correct snr value!)
   Eb_by_N0_bpsk_db(i) = 10*log10(Eb_by_N0_bpsk(i)); % convert to dB
   BER_bpsk(i) = 1./2.*erfc(sqrt(Eb_by_N0_bpsk(i)));
end
%%
figure(10)
plot(r,Eb_by_N0_bpsk_db, 'Linewidth', 3.5); 
title('Normalized SNR')
xlabel('distance(m)')
ylabel('Normalized SNR for BPSK')

%%
figure(11)
semilogy(r,BER_bpsk, 'Linewidth', 3.5)
title('BER BPSK')
xlabel('distance(m)')
ylabel('BER')
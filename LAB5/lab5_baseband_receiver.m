
load('myfile.mat')
load('OFDM.mat')
%%
OFDM_received = OFDM % you may not need this based on what you have loaded!

%% S/P converter
y_paral=reshape(OFDM_received,size(x_modulating_with_cpe,1),size(x_modulating_with_cpe,2)); % converting serial to parallel data

y_parall=conj(y_paral'); % reshaping rows and columns for processing
%% Remove cyclic prefix
y_recieved_without_cpe = y_parall(:,cps+1:end); % the cyclic prefix is at the start now...
%% FFT Block
for ii = 1:size(y_recieved_without_cpe,1)  
    y_fft(ii,:) = fft(y_recieved_without_cpe(ii,:)); 
end 
%% Remove flipped conj and padding
y_removed_fd = y_fft(:,2:size(x_par_pad,2)); % cleaning up our signal and getting rid of the 0 Hz carrier just as we had added only 0s to it in transmission
%% P/S converter
y_fd=conj(y_removed_fd'); % now back to proper shape of rows and columns
recvd_serial_data = reshape(y_fd,1,size(y_fd,1)*size(y_fd,2)); % this should be one row, 'X' columns data...


%% BPSK demodulation 
% we simply check which column of the symbol book the received signal
% stream corresponds to...
symbol_book_check = real(symbol_book');
recvd_serial_data_check =real(recvd_serial_data)';
rec_syms_col = knnsearch([symbol_book_check], [recvd_serial_data_check]); % this tells us the index of the column of symbol_book to which the recvd_signal_data belongs to 
% (column 1 = 0, column 2 = 1)
%%
LL=length(rec_syms_col);
Length=LL;

for kk=1:LL
if rec_syms_col(kk) == 1
    rec_syms(kk) = 0;
end
if rec_syms_col(kk) == 2
    rec_syms(kk) = 1;
end
end
%% error calculation % 
[Err,ratio] = biterr(stream,rec_syms)%find errors...
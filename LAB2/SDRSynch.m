%% Secction 1: Calibrates your SDR as a reciever
sampleRate = 200e3;
centerFreq = 2.42e9;
fRef = 80e3;
s1 = exp(1j*2*pi*20e3*[0:10000-1]'/sampleRate);
s2 = exp(1j*2*pi*40e3*[0:10000-1]'/sampleRate);
s3 = exp(1j*2*pi*fRef*[0:10000-1]'/sampleRate);
s = s1 + s2 + s3;
s = 0.6*s/max(abs(s));
numSamples = 1024*1024;
rx = sdrrx('Pluto', 'RadioID', 'usb:0', 'CenterFrequency', centerFreq, ...
 'BasebandSampleRate', sampleRate, 'SamplesPerFrame', numSamples, ...
 'OutputDataType', 'double', 'ShowAdvancedProperties', true);
% The info method shows the actual values of various hardware-related
% properties
rxRadioInfo = info(rx)
%% Section 2: Sets up the display
disp(['Capture signal and observe the frequency offset' newline])
receivedSig = rx();
%% Section 3: Performs FFT Decimation
y = fftshift(abs(fft(receivedSig)));
[~, idx] = findpeaks(y,'MinPeakProminence',max(0.5*y));
fReceived = (max(idx)-numSamples/2-1)/numSamples*sampleRate;
% Plot the spectrum
%% Section 4: Identifies the frequency offset
rx.FrequencyCorrection = ((fReceived - fRef)/(fRef + centerFreq))*10.^6;
%rx.FrequencyCorrection = ((centerFreq + fRef)/(fReceived + centerFreq))
msg = sprintf(['Based on the tone detected at %.3f kHz, ' ...
 'FrequencyCorrection of the receiver should be set to %.4f'], ...
 fReceived/1000, rx.FrequencyCorrection);
disp(msg);
rxRadioInfo = info(rx)
%% Section 5: Verifies the frequency offset
disp(['Capture signal and verify frequency correction' newline])
for i = 1:20 % We capture 20 samples here, but only use 1. Why?
 receivedSig = rx();
end
% Find the tone that corresponds to the 80 kHz transmitted tone
% fReceived2 should be very close to 80 kHz
y = fftshift(abs(fft(receivedSig)));
[~,idx] = findpeaks(y,'MinPeakProminence',max(0.5*y));
fReceived2 = (max(idx)-numSamples/2-1)/numSamples*sampleRate;

sa = dsp.SpectrumAnalyzer('SampleRate', sampleRate, 'SpectralAverages', 4);
sa.Title = sprintf('Tone Expected at 80 kHz, Actually Received at %.3f kHz', ...
 fReceived2/1000);
receivedSig = reshape(receivedSig, [], 16); % Reshape into 16 columns
for i = 1:size(receivedSig, 2)
 sa(receivedSig(:,i));
end





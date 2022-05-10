%% Section 1 Clears command prompt

clc;
clear;

%% Section 2 Sets the radio parameters

deviceName = 'Pluto';
samplerate = 528e3;
fmStationFrequency = 94.9e6; % adjust this to select another FM station

%% Section 3 Sets the parameters that you need to receive

rx = sdrrx(deviceName,'BasebandSampleRate',samplerate,'CenterFrequency',fmStationFrequency,'OutputDataType','double');

%% Section 4 Sets the window to capture

capture(rx,5,'Seconds','Filename','FMRecording.bb');

%% Section 5 Stops the interface with the SDR

release(rx);

%% Section 6 Sets the reading sample rate and reads from the bb file
bbr = comm.BasebandFileReader('FMRecording.bb');
bbr.SamplesPerFrame = 4400;

%% Section 7 Demodulates the signal

fmbDemod = comm.FMBroadcastDemodulator('AudioSampleRate',48e3, ... 
    'SampleRate',bbr.Metadata.BasebandSampleRate,'PlaySound',true);
while ~isDone(bbr)
        fmbDemod(bbr());
end

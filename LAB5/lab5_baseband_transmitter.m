clc;
clear all;
close all;

%% modulator 
%% Generated bits

stream=[-1,1,-1,1,1,1,-1,1,1,-1,1,1,1,-1,1,-1,1,-1,-1,1,1,-1,1,1,1,-1,-1,1,-1,-1,1,1,1,1,-1,-1,-1,-1,1,1,-1,1,-1,1,-1,1,-1,1,1,-1,1,-1,-1,1,-1,1,-1,1,1,1,-1,1,1,-1,-1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,-1,1,-1,-1,-1,1,1,1,1,1,-1,-1,-1,-1,1,1,1,1,-1,1,1,-1,1,1,1,1,-1,-1,1,-1,1,1,1,-1,-1,-1,-1,1,1,1,-1,-1,1,-1,1,1,-1,1,1,1,1,1,1,1,1,-1,-1,1,1,-1,1,1,1,-1,1,1,1,1,-1,1,-1,-1,1,1,1,-1,1,-1,1,1,1,-1,-1,-1,1,1,1,1,1,-1,1,-1,-1,1,1,1,-1,1,1,-1,-1,-1,1,-1,1,1,-1,-1,1,1,1,1,-1,-1,1,1,1,-1,1,-1,-1,-1,1,-1,-1,-1,-1,1,1,-1,-1,1,1,1,1,-1,1,-1,-1,-1,-1,-1,-1,-1,1,-1,-1,1,1,1,1,1,1,1,-1,-1,-1,-1,1,-1,-1,1,-1,1,1,1,1,-1,-1,1,1,1,1,1,1,-1,-1,-1,1,1,1,1,-1,1,1,-1,1,1,-1,1,1,-1,-1,-1,1,1,-1,1,1,1,1,-1,-1,1,-1,-1,-1,1,-1,1,-1,1,1,-1,1,-1,1,1,-1,-1,-1,-1,-1,1,1,1,1,-1,1,-1,1,-1,1,-1,1,1,1,1,1,-1,-1,1,1,1,-1,1,1,-1,1,1,-1,1,1,-1,1,1,1,1,-1,-1,1,1,-1,-1,-1,-1,-1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,-1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,-1,-1,-1,1,-1,-1,-1,-1,-1,1,-1,1,1,-1,1,1,1,1,-1,1,-1,-1,1,1,-1,-1,1,1,1,1,-1,1,-1,1,-1,1,1,-1,-1,-1,-1,-1,1,-1,1,-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,-1,1,1,1,-1,-1,1,1,1,1,1,1,-1,1,-1,-1,-1,-1,1,-1,1,-1,-1,1,-1,1,1,-1,-1,1,1,1,-1,-1,1,1,1,1,1,-1,1,-1,-1,1,-1,1,1,-1,-1,-1,1,-1,-1,-1,-1,-1,-1,-1,-1,1,-1,1,1,-1,-1,1,-1,-1,1,-1,1,1,1,1,-1,1,-1,1,1,1,-1,-1,-1,1,1,1,1,-1,1,1,-1,1,-1,-1,1,-1,1,1,1,1,1,1,-1,-1,1,-1,-1,-1,1,1,1,-1,-1,1,1,1,-1,-1,1,1,1,-1,1,-1,1,-1,-1,-1,1,1,-1,-1,-1,1,1,-1,-1,-1,1,1,-1,-1,1,-1,-1,1,-1,1,1,1,-1,1,1,-1,1,-1,-1,-1,1,1,-1,1,-1,1,1,-1,1,-1,1,-1,-1,1,1,-1,-1,1,-1,1,1,1,-1,-1,1,1,1,1,1,-1,1,-1,-1,-1,1,1,1,1,1,1,1,-1,-1,-1,1,1,1,1,-1,1,1,-1,1,-1,-1,-1,-1,1,-1,-1,-1,-1,-1,1,-1,-1,-1,1,-1,1,1,-1,1,1,-1,1,1,-1,-1,1,1,1,-1,1,1,1,1,-1,-1,-1,-1,-1,-1,1,1,-1,-1,-1,-1,-1,-1,1,-1,1,-1,-1,-1,-1,-1,1,-1,1,1,1,-1,-1,-1,1,1,-1,-1,1,1,-1,1,1,-1,1,-1,-1,1,1,-1,1,1,1,-1,-1,-1,1,-1,-1,1,1,-1,-1,-1,1,-1,-1,1,-1,1,1,-1,-1,1,1,1,1,-1,1,-1,-1,-1,1,1,1,1,-1,1,-1,1,1,1,-1,-1,1,1,1,1,1,1,1,-1,1,-1,-1,1,-1,-1,1,-1,-1,-1,-1,-1,1,-1,-1,-1,-1,1,-1,1,1,1,1,-1,-1,1,-1,1,1,1,-1,1,-1,-1,-1,1,1,1,1,-1,1,-1,1,1,1,1,1,-1,-1,1,-1,-1,-1,-1,1,1,1,-1,1,1,1,1,-1,1,-1,1,1,1,-1,1,-1,-1,1,1,-1,1,1,-1,-1,-1,-1,-1,1,-1,1,-1,1,-1,1,1,-1,-1,-1,-1,-1,-1,1,1,1,-1,1,-1,1,1,-1,1,1,1,1,1,-1,-1,-1,1,1,-1,-1,1,1,1,-1,1,1,-1,-1,-1,1,1,-1,1,1,-1,-1,1,1,1,1,-1,-1,1,1,1,-1,1,1,1,-1,1,1,-1,1,1,-1,1,-1,1,-1,-1,1,-1,-1,-1,-1,1,1,-1,-1,1,-1,1,1,1,1,-1,1,1,1,1,1,-1,1,1,-1,-1,1,-1,1,1,1,1,-1,-1,1,1,1,1,1,1,-1,1,1,-1,1,1,-1,1,-1,-1,-1,-1,-1,1,-1,-1,-1,1,-1,1,1,1,-1,1,1,-1,-1,-1,1,-1,1,1,1,-1,-1,1,-1,-1,-1,-1,1,1,1,1,-1,-1,-1,1,1,-1,-1,-1,1,1,-1,-1,-1,1,-1,1,-1,1,1,1,1,1,-1,1,-1,-1,-1,1,-1,1,1,1,1,-1,1,1,1,-1,1,1,1,-1,-1,-1,1,1,1,-1,1,1,1,-1,-1,1,1,1,-1,-1,1,1,-1,-1,-1,1,1,1,1,-1,1,1,-1,-1,-1,-1,1,1,1,-1,1,-1,1,1,1,-1,1,1,-1,1,1,1,-1,1,-1,1,-1,-1,1,1,-1,1,1,1,-1,-1,1,-1,-1,1,1,1,1,-1,-1,-1,-1,1,1,-1,1,-1,1,-1,1,-1,1,1,-1,1,-1,-1,1,-1,1,-1,1,1,1,-1,1,1,-1,-1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,-1,1,-1,-1,-1,1,1,1,1,1,-1,-1,-1,-1,1,1,1,1,-1,1,1,-1,1,1,1,1,-1,-1,1,-1,1,1,1,-1,-1,-1,-1,1,1,1,-1,-1,1,-1,1,1,-1,1,1,1,1,1,1,1,1,-1,-1,1,1,-1,1,1,1,-1,1,1,1,1,-1,1,-1,-1,1,1,1,-1,1,-1,1,1,1,-1,-1,-1,1,1,1,1,1,-1,1,-1,-1,1,1,1,-1,1,1,-1,-1,-1,1,-1,1,1,-1,-1,1,1,1,1,-1,-1,1,1,1,-1,1,-1,-1,-1,1,-1,-1,-1,-1,1,1,-1,-1,1,1,1,1,-1,1,-1,-1,-1,-1,-1,-1,-1,1,-1,-1,1,1,1,1,1,1,1,-1,-1,-1,-1,1,-1,-1,1,-1,1,1,1,1,-1,-1,1,1,1,1,1,1,-1,-1,-1,1,1,1,1,-1,1,1,-1,1,1,-1,1,1,-1,-1,-1,1,1,-1,1,1,1,1,-1,-1,1,-1,-1,-1,1,-1,1,-1,1,1,-1,1,-1,1,1,-1,-1,-1,-1,-1,1,1,1,1,-1,1,-1,1,-1,1,-1,1,1,1,1,1,-1,-1,1,1,1,-1,1,1,-1,1,1,-1,1,1,-1,1,1,1,1,-1,-1,1,1,-1,-1,-1,-1,-1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,-1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,-1,-1,-1,1,-1,-1,-1,-1,-1,1,-1,1,1,-1,1,1,1,1,-1,1,-1,-1,1,1,-1,-1,1,1,1,1,-1,1,-1,1,-1,1,1,-1,-1,-1,-1,-1,1,-1,1,-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,-1,1,1,1,-1,-1,1,1,1,1,1,1,-1,1,-1,-1,-1,-1,1,-1,1,-1,-1,1,-1,1,1,-1,-1,1,1,1,-1,-1,1,1,1,1,1,-1,1,-1,-1,1,-1,1,1,-1,-1,-1,1,-1,-1,-1,-1,-1,-1,-1,-1,1,-1,1,1,-1,-1,1,-1,-1,1,-1,1,1,1,1,-1,1,-1,1,1,1,-1,-1,-1,1,1,1,1,-1,1,1,-1,1,-1,-1,1,-1,1,1,1,1,1,1,-1,-1,1,-1,-1,-1,1,1,1,-1,-1,1,1,1,-1,-1,1,1,1,-1,1,-1,1,-1,-1,-1,1,1,-1,-1,-1,1,1,-1,-1,-1,1,1,-1,-1,1,-1,-1,1,-1,1,1,1,-1,1,1,-1,1,-1,-1,-1,1,1,-1,1,-1,1,1,-1,1,-1,1,-1,-1,1,1,-1,-1,1,-1,1,1,1,-1,-1,1,1,1,1,1,-1,1,-1,-1,-1,1,1,1,1,1,1,1,-1,-1,-1,1,1,1,1,-1,1,1,-1,1,-1,-1,-1,-1,1,-1,-1,-1,-1,-1,1,-1,-1,-1,1,-1,1,1,-1,1,1,-1,1,1,-1,-1,1,1,1,-1,1,1,1,1,-1,-1,-1,-1,-1,-1,1,1,-1,-1,-1,-1,-1,-1,1,-1,1,-1,-1,-1,-1,-1,1,-1,1,1,1,-1,-1,-1,1,1,-1,-1,1,1,-1,1,1,-1,1,-1,-1,1,1,-1,1,1,1,-1,-1,-1,1,-1,-1,1,1,-1,-1,-1,1,-1,-1,1,-1,1,1,-1,-1,1,1,1,1,-1,1,-1,-1,-1,1,1,1,1,-1,1,-1,1,1,1,-1,-1,1,1,1,1,1,1,1,-1,1,-1,-1,1,-1,-1,1,-1,-1,-1,-1,-1,1,-1,-1,-1,-1,1,-1,1,1,1,1,-1,-1,1,-1,1,1,1,-1,1,-1,-1,-1,1,1,1,1,-1,1,-1,1,1,1,1,1,-1,-1,1,-1,-1,-1,-1,1,1,1,-1,1,1,1,1,-1,1,-1,1,1,1,-1,1,-1,-1,1,1,-1,1,1,-1,-1,-1,-1,-1,1,-1,1,-1,1,-1,1,1,-1,-1,-1,-1,-1,-1,1,1,1,-1,1,-1,1,1,-1,1,1,1,1,1,-1,-1,-1,1,1,-1,-1,1,1,1,-1,1,1,-1,-1,-1,1,1,-1,1,1,-1,-1,1,1,1,1,-1,-1,1,1,1,-1,1,1,1,-1,1,1,-1,1,1,-1,1,-1,1,-1,-1,1,-1,-1,-1,-1,1,1,-1,-1,1,-1,1,1,1,1,-1,1,1,1,1,1,-1,1,1,-1,-1,1,-1,1,1,1,1,-1,-1,1,1,1,1,1,1,-1,1,1,-1,1,1,-1,1,-1,-1,-1,-1,-1,1,-1,-1,-1,1,-1,1,1,1,-1,1,1,-1,-1,-1,1,-1,1,1,1,-1,-1,1,-1,-1,-1,-1,1,1,1,1,-1,-1,-1,1,1,-1,-1,-1,1,1,-1,-1,-1,1,-1,1,-1,1,1,1,1,1,-1,1,-1,-1,-1,1,-1,1,1,1,1,-1,1,1,1,-1,1,1,1,-1,-1,-1,1,1,1,-1,1,1,1,-1,-1,1,1,1,-1,-1,1,1,-1,-1,-1,1,1,1,1,-1,1,1,-1,-1,-1,-1,1,1,1,-1,1,-1,1,1,1,-1,1,1,-1,1,1,1,-1,1,-1,1,-1,-1,1,1,-1,1,1,1,-1,-1,1,-1,-1,1,1,1,1,-1,-1,-1,-1,1,1,-1,1,-1,1,-1,1,-1,1,1,-1,1,-1,-1,1,-1,1,-1,1,1,1,-1,1,1,-1,-1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,-1,1,-1,-1,-1,1,1,1,1,1,-1,-1,-1,-1,1,1,1,1,-1,1,1,-1,1,1,1,1,-1,-1,1,-1,1,1,1,-1,-1,-1,-1,1,1,1,-1,-1,1,-1,1,1,-1,1,1,1,1,1,1,1,1,-1,-1,1,1,-1,1,1,1,-1,1,1,1,1,-1,1,-1,-1,1,1,1,-1,1,-1,1,1,1,-1,-1,-1,1,1,1,1,1,-1,1,-1,-1,1,1,1,-1,1,1,-1,-1,-1,1,-1,1,1,-1,-1,1,1,1,1,-1,-1,1,1,1,-1,1,-1,-1,-1,1,-1,-1,-1,-1,1,1,-1,-1,1,1,1,1,-1,1,-1,-1,-1,-1,-1,-1,-1,1,-1,-1,1,1,1,1,1,1,1,-1,-1,-1,-1,1,-1,-1,1,-1,1,1,1,1,-1,-1,1,1,1,1,1,1,-1,-1,-1,1,1,1,1,-1,1,1,-1,1,1,-1,1,1,-1,-1,-1,1,1,-1,1,1,1,1,-1,-1,1,-1,-1,-1,1,-1,1,-1,1,1,-1,1,-1,1,1,-1,-1,-1,-1,-1,1,1,1,1,-1,1,-1,1,-1,1,-1,1,1,1,1,1,-1,-1,1,1,1,-1,1,1,-1,1,1,-1,1,1,-1,1,1,1,1,-1,-1,1,1,-1,-1,-1,-1,-1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,-1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,-1,-1,-1,1,-1,-1,-1,-1,-1,1,-1,1,1,-1,1,1,1,1,-1,1,-1,-1,1,1,-1,-1,1,1,1,1,-1,1,-1,1,-1,1,1,-1,-1,-1,-1,-1,1,-1,1,-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,-1,1,1,1,-1,-1,1,1,1,1,1,1,-1,1,-1,-1,-1,-1,1,-1,1,-1,-1,1,-1,1,1,-1,-1,1,1,1,-1,-1,1,1,1,1,1,-1,1,-1,-1,1,-1,1,1,-1,-1,-1,1,-1,-1,-1,-1,-1,-1,-1,-1,1,-1,1,1,-1,-1,1,-1,-1,1,-1,1,1,1,1,-1,1,-1,1,1,1,-1,-1,-1,1,1,1,1,-1,1,1,-1,1,-1,-1,1,-1,1,1,1,1,1,1,-1,-1,1,-1,-1,-1,1,1,1,-1,-1,1,1,1,-1,-1,1,1,1,-1,1,-1,1,-1,-1,-1,1,1,-1,-1,-1,1,1,-1,-1,-1,1,1,-1,-1,1,-1,-1,1,-1,1,1,1,-1,1,1,-1,1,-1,-1,-1,1,1,-1,1,-1,1,1,-1,1,-1,1,-1,-1,1,1,-1,-1,1,-1,1,1,1,-1,-1,1,1,1,1,1,-1,1,-1,-1,-1,1,1,1,1,1,1,1,-1,-1,-1,1,1,1,1,-1,1,1,-1,1,-1,-1,-1,-1,1,-1,-1,-1,-1,-1,1,-1,-1,-1,1,-1,1,1,-1,1,1,-1,1,1,-1,-1,1,1,1,-1,1,1,1,1,-1,-1,-1,-1,-1,-1,1,1,-1,-1,-1,-1,-1,-1,1,-1,1,-1,-1,-1,-1,-1,1,-1,1,1,1,-1,-1,-1,1,1,-1,-1,1,1,-1,1,1,-1,1,-1,-1,1,1,-1,1,1,1,-1,-1,-1,1,-1,-1,1,1,-1,-1,-1,1,-1,-1,1,-1,1,1,-1,-1,1,1,1,1,-1,1,-1,-1,-1,1,1,1,1,-1,1,-1,1,1,1,-1,-1,1,1,1,1,1,1,1,-1,1,-1,-1,1,-1,-1,1,-1,-1,-1,-1,-1,1,-1,-1,-1,-1,1,-1,1,1,1,1,-1,-1,1,-1,1,1,1,-1,1,-1,-1,-1,1,1,1,1,-1,1,-1,1,1,1,1,1,-1,-1,1,-1,-1,-1,-1,1,1,1,-1,1,1,1,1,-1,1,-1,1,1,1,-1,1,-1,-1,1,1,-1,1,1,-1,-1,-1,-1,-1,1,-1,1,-1,1,-1,1,1,-1,-1,-1,-1,-1,-1,1,1,1,-1,1,-1,1,1,-1,1,1,1,1,1,-1,-1,-1,1,1,-1,-1,1,1,1,-1,1,1,-1,-1,-1,1,1,-1,1,1,-1,-1,1,1,1,1,-1,-1,1,1,1,-1,1,1,1,-1,1,1,-1,1,1,-1,1,-1,1,-1,-1,1,-1,-1,-1,-1,1,1,-1,-1,1,-1,1,1,1,1,-1,1,1,1,1,1,-1,1,1,-1,-1,1,-1,1,1,1,1,-1,-1,1,1,1,1,1,1,-1,1,1,-1,1,1,-1,1,-1,-1,-1,-1,-1,1,-1,-1,-1,1,-1,1,1,1,-1,1,1,-1,-1,-1,1,-1,1,1,1,-1,-1,1,-1,-1,-1,-1,1,1,1,1,-1,-1,-1,1,1,-1,-1,-1,1,1,-1,-1,-1,1,-1,1,-1,1,1,1,1,1,-1,1,-1,-1,-1,1,-1,1,1,1,1,-1,1,1,1,-1,1,1,1,-1,-1,-1,1,1,1,-1,1,1,1,-1,-1,1,1,1,-1,-1,1,1,-1,-1,-1,1,1,1,1,-1,1,1,-1,-1,-1,-1,1,1,1,-1,1,-1,1,1,1,-1,1,1,-1,1,1,1,-1,1,-1,1,-1,-1,1,1,-1,1,1,1,-1,-1,1,-1,-1,1,1,1,1,-1,-1,-1,-1,1,1,-1,1,-1,1,-1,1,-1,1,1,-1,1,-1,-1,1,-1,1,-1,1,1,1,-1,1,1,-1,-1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,-1,1,-1,-1,-1,1,1,1,1,1,-1,-1,-1,-1,1,1,1,1,-1,1,1,-1,1,1,1,1,-1,-1,1,-1,1,1,1,-1,-1,-1,-1,1,1,1,-1,-1,1,-1,1,1,-1,1,1,1,1,1,1,1,1,-1,-1,1,1,-1,1,1,1,-1,1,1,1,1,-1,1,-1,-1,1,1,1,-1,1,-1,1,1,1,-1,-1,-1,1,1,1,1,1,-1,1,-1,-1,1,1,1,-1,1,1,-1,-1,-1,1,-1,1,1,-1,-1,1,1,1,1,-1,-1,1,1,1,-1,1,-1,-1,-1,1,-1,-1,-1,-1,1,1,-1,-1,1,1,1,1,-1,1,-1,-1,-1,-1,-1,-1,-1,1,-1,-1,1,1,1,1,1,1,1,-1,-1,-1,-1,1,-1,-1,1,-1,1,1,1,1,-1,-1,1,1,1,1,1,1,-1,-1,-1,1,1,1,1,-1,1,1,-1,1,1,-1,1,1,-1,-1,-1,1,1,-1,1,1,1,1,-1,-1,1,-1,-1,-1,1,-1,1,-1,1,1,-1,1,-1,1,1,-1,-1,-1,-1,-1,1,1,1,1,-1,1,-1,1,-1,1,-1,1,1,1,1,1,-1,-1,1,1,1,-1,1,1,-1,1,1,-1,1,1,-1,1,1,1,1,-1,-1,1,1,-1,-1,-1,-1,-1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,-1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,-1,-1,-1,1,-1,-1,-1,-1,-1,1,-1,1,1,-1,1,1,1,1,-1,1,-1,-1,1,1,-1,-1,1,1,1,1,-1,1,-1,1,-1,1,1,-1,-1,-1,-1,-1,1,-1,1,-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,-1,1,1,1,-1,-1,1,1,1,1,1,1,-1,1,-1,-1,-1,-1,1,-1,1,-1,-1,1,-1,1,1,-1,-1,1,1,1,-1,-1,1,1,1,1,1,-1,1,-1,-1,1,-1,1,1,-1,-1,-1,1,-1,-1,-1,-1,-1,-1,-1,-1,1,-1,1,1,-1,-1,1,-1,-1,1,-1,1,1,1,1,-1,1,-1,1,1,1,-1,-1,-1,1,1,1,1,-1,1,1,-1,1,-1,-1,1,-1,1,1,1,1,1,1,-1,-1,1,-1,-1,-1,1,1,1,-1,-1,1,1,1,-1,-1,1,1,1,-1,1,-1,1,-1,-1,-1,1,1,-1,-1,-1,1,1,-1,-1,-1,1,1,-1,-1,1,-1,-1,1,-1,1,1,1,-1,1,1,-1,1,-1,-1,-1,1,1,-1,1,-1,1,1,-1,1,-1,1,-1,-1,1,1,-1,-1,1,-1,1,1,1,-1,-1,1,1,1,1,1,-1,1,-1,-1,-1,1,1,1,1,1,1,1,-1,-1,-1,1,1,1,1,-1,1,1,-1,1,-1,-1,-1,-1,1,-1,-1,-1,-1,-1,1,-1,-1,-1,1,-1,1,1,-1,1,1,-1,1,1,-1,-1,1,1,1,-1,1,1,1,1,-1,-1,-1,-1,-1,-1,1,1,-1,-1,-1,-1,-1,-1,1,-1,1,-1,-1,-1,-1,-1,1,-1,1,1,1,-1,-1,-1,1,1,-1,-1,1,1,-1,1,1,-1,1,-1,-1,1,1,-1,1,1,1,-1,-1,-1,1,-1,-1,1,1,-1,-1,-1,1,-1,-1,1,-1,1,1,-1,-1,1,1,1,1,-1,1,-1,-1,-1,1,1,1,1,-1,1,-1,1,1,1,-1,-1,1,1,1,1,1,1,1,-1,1,-1,-1,1,-1,-1,1,-1,-1,-1,-1,-1,1,-1,-1,-1,-1,1,-1,1,1,1,1,-1,-1,1,-1,1,1,1,-1,1,-1,-1,-1,1,1,1,1,-1,1,-1,1,1,1,1,1,-1,-1,1,-1,-1,-1,-1,1,1,1,-1,1,1,1,1,-1,1,-1,1,1,1,-1,1,-1,-1,1,1,-1,1,1,-1,-1,-1,-1,-1,1,-1,1,-1,1,-1,1,1,-1,-1,-1,-1,-1,-1,1,1,1,-1,1,-1,1,1,-1,1,1,1,1,1,-1,-1,-1,1,1,-1,-1,1,1,1,-1,1,1,-1,-1,-1,1,1,-1,1,1,-1,-1,1,1,1,1,-1,-1,1,1,1,-1,1,1,1,-1,1,1,-1,1,1,-1,1,-1,1,-1,-1,1,-1,-1,-1,-1,1,1,-1,-1,1,-1,1,1,1,1,-1,1,1,1,1,1,-1,1,1,-1,-1,1,-1,1,1,1,1,-1,-1,1,1,1,1,1,1,-1,1,1,-1,1,1,-1,1,-1,-1,-1,-1,-1,1,-1,-1,-1,1,-1,1,1,1,-1,1,1,-1,-1,-1,1,-1,1,1,1,-1,-1,1,-1,-1,-1,-1,1,1,1,1,-1,-1,-1,1,1,-1,-1,-1,1,1,-1,-1,-1,1,-1,1,-1,1,1,1,1,1,-1,1,-1,-1,-1,1,-1,1,1,1,1,-1,1,1,1,-1,1,1,1,-1,-1,-1,1,1,1,-1,1,1,1,-1,-1,1,1,1,-1,-1,1,1,-1,-1,-1,1,1,1,1,-1,1,1,-1,-1,-1,-1,1,1,1,-1,1,-1,1,1,1,-1,1,1,-1,1,1,1,-1,1,-1,1,-1,-1,1,1,-1,1,1,1,-1,-1,1,-1,-1,1,1,1,1,-1,-1,-1,-1,1,1,-1,1,-1,1,-1,1,-1,1,1,-1,1,-1,-1,1,-1,1,-1,1,1,1,-1,1,1,-1,-1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,-1,1,-1,-1,-1,1,1,1,1,1,-1,-1,-1,-1,1,1,1,1,-1,1,1,-1,1,1,1,1,-1,-1,1,-1,1,1,1,-1,-1,-1,-1,1,1,1,-1,-1,1,-1,1,1,-1,1,1,1,1,1,1,1,1,-1,-1,1,1,-1,1,1,1,-1,1,1,1,1,-1,1,-1,-1,1,1,1,-1,1,-1,1,1,1,-1,-1,-1,1,1,1,1,1,-1,1,-1,-1,1,1,1,-1,1,1,-1,-1,-1,1,-1,1,1,-1,-1,1,1,1,1,-1,-1,1,1,1,-1,1,-1,-1,-1,1,-1,-1,-1,-1,1,1,-1,-1,1,1,1,1,-1,1,-1,-1,-1,-1,-1,-1,-1,1,-1,-1,1,1,1,1,1,1,1,-1,-1,-1,-1,1,-1,-1,1,-1,1,1,1,1,-1,-1,1,1,1,1,1,1,-1,-1,-1,1,1,1,1,-1,1,1,-1,1,1,-1,1,1,-1,-1,-1,1,1,-1,1,1,1,1,-1,-1,1,-1,-1,-1,1,-1,1,-1,1,1,-1,1,-1,1,1,-1,-1,-1,-1,-1,1,1,1,1,-1,1,-1,1,-1,1,-1,1,1,1,1,1,-1,-1,1,1,1,-1,1,1,-1,1,1,-1,1,1,-1,1,1,1,1,-1,-1,1,1,-1,-1,-1,-1,-1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,-1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,-1,-1,-1,1,-1,-1,-1,-1,-1,1,-1,1,1,-1,1,1,1,1,-1,1,-1,-1,1,1,-1,-1,1,1,1,1,-1,1,-1,1,-1,1,1,-1,-1,-1,-1,-1,1,-1,1,-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,-1,1,1,1,-1,-1,1,1,1,1,1,1,-1,1,-1,-1,-1,-1,1,-1,1,-1,-1,1,-1,1,1,-1,-1,1,1,1,-1,-1,1,1,1,1,1,-1,1,-1,-1,1,-1,1,1,-1,-1,-1,1,-1,-1,-1,-1,-1,-1,-1,-1,1,-1,1,1,-1,-1,1,-1,-1,1,-1,1,1,1,1,-1,1,-1,1,1,1,-1,-1,-1,1,1,1,1,-1,1,1,-1,1,-1,-1,1,-1,1,1,1,1,1,1,-1,-1,1,-1,-1,-1,1,1,1,-1,-1,1,1,1,-1,-1,1,1,1,-1,1,-1,1,-1,-1,-1,1,1,-1,-1,-1,1,1,-1,-1,-1,1,1,-1,-1,1,-1,-1,1,-1,1,1,1,-1,1,1,-1,1,-1,-1,-1,1,1,-1,1,-1,1,1,-1,1,-1,1,-1,-1,1,1,-1,-1,1,-1,1,1,1,-1,-1,1,1,1,1,1,-1,1,-1,-1,-1,1,1,1,1,1,1,1,-1,-1,-1,1,1,1,1,-1,1,1,-1,1,-1,-1,-1,-1,1,-1,-1,-1,-1,-1,1,-1,-1,-1,1,-1,1,1,-1,1,1,-1,1,1,-1,-1,1,1,1,-1,1,1,1,1,-1,-1,-1,-1,-1,-1,1,1,-1,-1,-1,-1,-1,-1,1,-1,1,-1,-1,-1,-1,-1,1,-1,1,1,1,-1,-1,-1,1,1,-1,-1,1,1,-1,1,1,-1,1,-1,-1,1,1,-1,1,1,1,-1,-1,-1,1,-1,-1,1,1,-1,-1,-1,1,-1,-1,1,-1,1,1,-1,-1,1,1,1,1,-1,1,-1,-1,-1,1,1,1,1,-1,1,-1,1,1,1,-1,-1,1,1,1,1,1,1,1,-1,1,-1,-1,1,-1,-1,1,-1,-1,-1,-1,-1,1,-1,-1,-1,-1,1,-1,1,1,1,1,-1,-1,1,-1,1,1,1,-1,1,-1,-1,-1,1,1,1,1,-1,1,-1,1,1,1,1,1,-1,-1,1,-1,-1,-1,-1,1,1,1,-1,1,1,1,1,-1,1,-1,1,1,1,-1,1,-1,-1,1,1,-1,1,1,-1,-1,-1,-1,-1,1,-1,1,-1,1,-1,1,1,-1,-1,-1,-1,-1,-1,1,1,1,-1,1,-1,1,1,-1,1,1,1,1,1,-1,-1,-1,1,1,-1,-1,1,1,1,-1,1,1,-1,-1,-1,1,1,-1,1,1,-1,-1,1,1,1,1,-1,-1,1,1,1,-1,1,1,1,-1,1,1,-1,1,1,-1,1,-1,1,-1,-1,1,-1,-1,-1,-1,1,1,-1,-1,1,-1,1,1,1,1,-1,1,1,1,1,1,-1,1,1,-1,-1,1,-1,1,1,1,1,-1,-1,1,1,1,1,1,1,-1,1,1,-1,1,1,-1,1,-1,-1,-1,-1,-1,1,-1,-1,-1,1,-1,1,1,1,-1,1,1,-1,-1,-1,1,-1,1,1,1,-1,-1,1,-1,-1,-1,-1,1,1,1,1,-1,-1,-1,1,1,-1,-1,-1,1,1,-1,-1,-1,1,-1,1,-1,1,1,1,1,1,-1,1,-1,-1,-1,1,-1,1,1,1,1,-1,1,1,1,-1,1,1,1,-1,-1,-1,1,1,1,-1,1,1,1,-1,-1,1,1,1,-1,-1,1,1,-1,-1,-1,1,1,1,1,-1,1,1,-1,-1,-1,-1,1,1,1,-1,1,-1,1,1,1,-1,1,1,-1,1,1,1,-1,1,-1,1,-1,-1,1,1,-1,1,1,1,-1,-1,1,-1,-1,1,1,1,1,-1,-1,-1,-1,1,1,-1,1,-1,1,-1,1,-1,1,1,-1,1,-1,-1,1,-1,1,-1,1,1,1,-1,1,1,-1,-1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,-1,1,-1,-1,-1,1,1,1,1,1,-1,-1,-1,-1,1,1,1,1,-1,1,1,-1,1,1,1,1,-1,-1,1,-1,1,1,1,-1,-1,-1,-1,1,1,1,-1,-1,1,-1,1,1,-1,1,1,1,1,1,1,1,1,-1,-1,1,1,-1,1,1,1,-1,1,1,1,1,-1,1,-1,-1,1,1,1,-1,1,-1,1,1,1,-1,-1,-1,1,1,1,1,1,-1,1,-1,-1,1,1,1,-1,1,1,-1,-1,-1,1,-1,1,1,-1,-1,1,1,1,1,-1,-1,1,1,1,-1,1,-1,-1,-1,1,-1,-1,-1,-1,1,1,-1,-1,1,1,1,1,-1,1,-1,-1,-1,-1,-1,-1,-1,1,-1,-1,1,1,1,1,1,1,1,-1,-1,-1,-1,1,-1,-1,1,-1,1,1,1,1,-1,-1,1,1,1,1,1,1,-1,-1,-1,1,1,1,1,-1,1,1,-1,1,1,-1,1,1,-1,-1,-1,1,1,-1,1,1,1,1,-1,-1,1,-1,-1,-1,1,-1,1,-1,1,1,-1,1,-1,1,1,-1,-1,-1,-1,-1,1,1,1,1,-1,1,-1,1,-1,1,-1,1,1,1,1,1,-1,-1,1,1,1,-1,1,1,-1,1,1,-1,1,1,-1,1,1,1,1,-1,-1,1,1,-1,-1,-1,-1,-1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,-1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,-1,-1,-1,1,-1,-1,-1,-1,-1,1,-1,1,1,-1,1,1,1,1,-1,1,-1,-1,1,1,-1,-1,1,1,1,1,-1,1,-1,1,-1,1,1,-1,-1,-1,-1,-1,1,-1,1,-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,-1,1,1,1,-1,-1,1,1,1,1,1,1,-1,1,-1,-1,-1,-1,1,-1,1,-1,-1,1,-1,1,1,-1,-1,1,1,1,-1,-1,1,1,1,1,1,-1,1,-1,-1,1,-1,1,1,-1,-1,-1,1,-1,-1,-1,-1,-1,-1,-1,-1,1,-1,1,1,-1,-1,1,-1,-1,1,-1,1,1,1,1,-1,1,-1,1,1,1,-1,-1,-1,1,1,1,1,-1,1,1,-1,1,-1,-1,1,-1,1,1,1,1,1,1,-1,-1,1,-1,-1,-1,1,1,1,-1,-1,1,1,1,-1,-1,1,1,1,-1,1,-1,1,-1,-1,-1,1,1,-1,-1,-1,1,1,-1,-1,-1,1,1,-1,-1,1,-1,-1,1,-1,1,1,1,-1,1,1,-1,1,-1,-1,-1,1,1,-1,1,-1,1,1,-1,1,-1,1,-1,-1,1,1,-1,-1,1,-1,1,1,1,-1,-1,1,1,1,1,1,-1,1,-1,-1,-1,1,1,1,1,1,1,1,-1,-1,-1,1,1,1,1,-1,1,1,-1,1,-1,-1,-1,-1,1,-1,-1,-1,-1,-1,1,-1,-1,-1,1,-1,1,1,-1,1,1,-1,1,1,-1,-1,1,1,1,-1,1,1,1,1,-1,-1,-1,-1,-1,-1,1,1,-1,-1,-1,-1,-1,-1,1,-1,1,-1,-1,-1,-1,-1,1,-1,1,1,1,-1,-1,-1,1,1,-1,-1,1,1,-1,1,1,-1,1,-1,-1,1,1,-1,1,1,1,-1,-1,-1,1,-1,-1,1,1,-1,-1,-1,1,-1,-1,1,-1,1,1,-1,-1,1,1,1,1,-1,1,-1,-1,-1,1,1,1,1,-1,1,-1,1,1,1,-1,-1,1,1,1,1,1,1,1,-1,1,-1,-1,1,-1,-1,1,-1,-1,-1,-1,-1,1,-1,-1,-1,-1,1,-1,1,1,1,1,-1,-1,1,-1,1,1,1,-1,1,-1,-1,-1,1,1,1,1,-1,1,-1,1,1,1,1,1,-1,-1,1,-1,-1,-1,-1,1,1,1,-1,1,1,1,1,-1,1,-1,1,1,1,-1,1,-1,-1,1,1,-1,1,1,-1,-1,-1,-1,-1,1,-1,1,-1,1,-1,1,1,-1,-1,-1,-1,-1,-1,1,1,1,-1,1,-1,1,1,-1,1,1,1,1,1,-1,-1,-1,1,1,-1,-1,1,1,1,-1,1,1,-1,-1,-1,1,1,-1,1,1,-1,-1,1,1,1,1,-1,-1,1,1,1,-1,1,1,1,-1,1,1,-1,1,1,-1,1,-1,1,-1,-1,1,-1,-1,-1,-1,1,1,-1,-1,1,-1,1,1,1,1,-1,1,1,1,1,1,-1,1,1,-1,-1,1,-1,1,1,1,1,-1,-1,1,1,1,1,1,1,-1,1,1,-1,1,1,-1,1,-1,-1,-1,-1,-1,1,-1,-1,-1,1,-1,1,1,1,-1,1,1,-1,-1,-1,1,-1,1,1,1,-1,-1,1,-1,-1,-1,-1,1,1,1,1,-1,-1,-1,1,1,-1,-1,-1,1,1,-1,-1,-1,1,-1,1,-1,1,1,1,1,1,-1,1,-1,-1,-1,1,-1,1,1,1,1,-1,1,1,1,-1,1,1,1,-1,-1,-1,1,1,1,-1,1,1,1,-1,-1,1,1,1,-1,-1,1,1,-1,-1,-1,1,1,1,1,-1,1,1,-1,-1,-1,-1,1,1,1,-1,1,-1,1,1,1,-1,1,1,-1,1,1,1,-1,1,-1,1,-1,-1,1,1,-1,1,1,1,-1,-1,1,-1,-1,1,1,1,1,-1,-1,-1,-1,1,1,-1,1,-1,1,-1,1,-1,1,1,-1,1,-1,-1,1,-1,1,-1,1,1,1,-1,1,1,-1,-1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,-1,1,-1,-1,-1,1,1,1,1,1,-1,-1,-1,-1,1,1,1,1,-1,1,1,-1,1,1,1,1,-1,-1,1,-1,1,1,1,-1,-1,-1,-1,1,1,1,-1,-1,1,-1,1,1,-1,1,1,1,1,1,1,1,1,-1,-1,1,1,-1,1,1,1,-1,1,1,1,1,-1,1,-1,-1,1,1,1,-1,1,-1,1,1,1,-1,-1,-1,1,1,1,1,1,-1,1,-1,-1,1,1,1,-1,1,1,-1,-1,-1,1,-1,1,1,-1,-1,1,1,1,1,-1,-1,1,1,1,-1,1,-1,-1,-1,1,-1,-1,-1,-1,1,1,-1,-1,1,1,1,1,-1,1,-1,-1,-1,-1,-1,-1,-1,1,-1,-1,1,1,1,1,1,1,1,-1,-1,-1,-1,1,-1,-1,1,-1,1,1,1,1,-1,-1,1,1,1,1,1,1,-1,-1,-1,1,1,1,1,-1,1,1,-1,1,1,-1,1,1,-1,-1,-1,1,1,-1,1,1,1,1,-1,-1,1,-1,-1,-1,1,-1,1,-1,1,1,-1,1,-1,1,1,-1,-1,-1,-1,-1,1,1,1,1,-1,1,-1,1,-1,1,-1,1,1,1,1,1,-1,-1,1,1,1,-1,1,1,-1,1,1,-1,1,1,-1,1,1,1,1,-1,-1,1,1,-1,-1,-1,-1,-1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,-1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,-1,-1,-1,1,-1,-1,-1,-1,-1,1,-1,1,1,-1,1,1,1,1,-1,1,-1,-1,1,1,-1,-1,1,1,1,1,-1,1,-1,1,-1,1,1,-1,-1,-1,-1,-1,1,-1,1,-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,-1,1,1,1,-1,-1,1,1,1,1,1,1,-1,1,-1,-1,-1,-1,1,-1,1,-1,-1,1,-1,1,1,-1,-1,1,1,1,-1,-1,1,1,1,1,1,-1,1,-1,-1,1,-1,1,1,-1,-1,-1,1,-1,-1,-1,-1,-1,-1,-1,-1,1,-1,1,1,-1,-1,1,-1,-1,1,-1,1,1,1,1,-1,1,-1,1,1,1,-1,-1,-1,1,1,1,1,-1,1,1,-1,1,-1,-1,1,-1,1,1,1,1,1,1,-1,-1,1,-1,-1,-1,1,1,1,-1,-1,1,1,1,-1,-1,1,1,1,-1,1,-1,1,-1,-1,-1,1,1,-1,-1,-1,1,1,-1,-1,-1,1,1,-1,-1,1,-1,-1,1,-1,1,1,1,-1,1,1,-1,1,-1,-1,-1,1,1,-1,1,-1,1,1,-1,1,-1,1,-1,-1,1,1,-1,-1,1,-1,1,1,1,-1,-1,1,1,1,1,1,-1,1,-1,-1,-1,1,1,1,1,1,1,1,-1,-1,-1,1,1,1,1,-1,1,1,-1,1,-1,-1,-1,-1,1,-1,-1,-1,-1,-1,1,-1,-1,-1,1,-1,1,1,-1,1,1,-1,1,1,-1,-1,1,1,1,-1,1,1,1,1,-1,-1,-1,-1,-1,-1,1,1,-1,-1,-1,-1,-1,-1,1,-1,1,-1,-1,-1,-1,-1,1,-1,1,1,1,-1,-1,-1,1,1,-1,-1,1,1,-1,1,1,-1,1,-1,-1,1,1,-1,1,1,1,-1,-1,-1,1,-1,-1,1,1,-1,-1,-1,1,-1,-1,1,-1,1,1,-1,-1,1,1,1,1,-1,1,-1,-1,-1,1,1,1,1,-1,1,-1,1,1,1,-1,-1,1,1,1,1,1,1,1,-1,1,-1,-1,1,-1,-1,1,-1,-1,-1,-1,-1,1,-1,-1,-1,-1,1,-1,1,1,1,1,-1,-1,1,-1,1,1,1,-1,1,-1,-1,-1,1,1,1,1,-1,1,-1,1,1,1,1,1,-1,-1,1,-1,-1,-1,-1,1,1,1,-1,1,1,1,1,-1,1,-1,1,1,1,-1,1,-1,-1,1,1,-1,1,1,-1,-1,-1,-1,-1,1,-1,1,-1,1,-1,1,1,-1,-1,-1,-1,-1,-1,1,1,1,-1,1,-1,1,1,-1,1,1,1,1,1,-1,-1,-1,1,1,-1,-1,1,1,1,-1,1,1,-1,-1,-1,1,1,-1,1,1,-1,-1,1,1,1,1,-1,-1,1,1,1,-1,1,1,1,-1,1,1,-1,1,1,-1,1,-1,1,-1,-1,1,-1,-1,-1,-1,1,1,-1,-1,1,-1,1,1,1,1,-1,1,1,1,1,1,-1,1,1,-1,-1,1,-1,1,1,1,1,-1,-1,1,1,1,1,1,1,-1,1,1,-1,1,1,-1,1,-1,-1,-1,-1,-1,1,-1,-1,-1,1,-1,1,1,1,-1,1,1,-1,-1,-1,1,-1,1,1,1,-1,-1,1,-1,-1,-1,-1,1,1,1,1,-1,-1,-1,1,1,-1,-1,-1,1,1,-1,-1,-1,1,-1,1,-1,1,1,1,1,1,-1,1,-1,-1,-1,1,-1,1,1,1,1,-1,1,1,1,-1,1,1,1,-1,-1,-1,1,1,1,-1,1,1,1,-1,-1,1,1,1,-1,-1,1,1,-1,-1,-1,1,1,1,1,-1,1,1,-1,-1,-1,-1,1,1,1,1,-1,-1,1,1,1];
stream=stream(1:4185);

LL=length(stream);
Length=LL;

for kk=1:LL
if stream(kk) == -1
    stream(kk) = 0;
end
end

%% Parameters
N_channels = 10; % Number of channels in OFDM 
cps=floor((1/8)*2*N_channels); % Cyclic prefix extension length  
fs=60e6; % [Hz] Sampling frequency in baseband
symbol_rate=2e6; % [Symbols/second] Frequency of the symbol or IF frequency 
M=2; %modulation level

bits_symbol = log2(M); % No of bits per symbol 

N_symbol=Length/bits_symbol; % This needs to be an integer, for simplicity 


ts = 1/fs; % Sampling time 


tc=1/symbol_rate; 


N_samples= fs/symbol_rate;% No of samples per symbol 


t_symbol=N_samples/fs; % 1 symbol duration


t_total=0:ts:N_symbol*N_samples*ts-ts; %total duration of transmitted symbols



%% BPSK modulation
mod_order = 1;
mod_ind = 2^(mod_order-1);
n = 0:pi/mod_ind:2*pi-pi/mod_ind;
in_phase = cos(n);
quadrature = 0;

symbol_book = (in_phase + quadrature*1i);
xn = symbol_book(stream+1); % goes through length of stream and makes equivalent for 0 and 1 


%% S/P converter

N_channels1=N_channels-1;

if mod(Length,N_channels1)==0
    xnn=xn;
else
    xnn=[xn zeros(1,(N_channels1-mod(Length,N_channels1)))] ; %padding
end

x_par=reshape(xnn, N_channels1, length(xnn)/N_channels1); %converting from serial to parallel
x_par=conj(x_par');
%% Padding and N-point flipped conjugate
x_par_pad = [zeros(size(x_par,1),1) x_par]; %adding zeros because?  
% for ii = 1:size(x_par_pad,1) 
%      x_par_pad_flipped(ii,:) = [x_par_pad(ii,:) 0 conj(flip(x_par_pad(ii,2:N_channels)))]; 
% end 
%% IFFT block
for ii = 1:size(x_par_pad,1)  
     IFFT_subcarriers(ii,:) = ifft(x_par_pad(ii,:)); %IFFT of each row is taken separately 
end 

%Plots
%%
figure(1)
L = length(t_total); 
NFFT = 2^nextpow2(L); 
f = fs*linspace(0,1,NFFT); 

X = fft(IFFT_subcarriers,NFFT)/L; 
S_X = 1./fs*L.*(X(1:NFFT).*X(1:NFFT)); 
S_X = (S_X); 
plot(f/1e9,S_X) 

title([num2str(length(stream)) ' symbols ' num2str(N_channels) ' channels']); 
xlabel('Frequency [GHz]') 
ylabel('S_X(f)') 
grid on 



%% Add cyclic prefix
x_modulating_with_cpe = [IFFT_subcarriers(:,end-cps+1:end), IFFT_subcarriers]; 
%% P/S converter
x_modulating_with_cpe=conj(x_modulating_with_cpe'); % contains all the information
OFDM=reshape(x_modulating_with_cpe, 1, size(x_modulating_with_cpe,1)*size(x_modulating_with_cpe,2)); % and we're done!

figure(2)
plot(OFDM);
%
t_new_2 = (0:1:(length(OFDM)-1))*ts;
%
figure(3)
plot(t_new_2,OFDM)
title('Modulated signal in Baseband'); 
xlabel('Time [s]') 
ylabel('Amplitude') 
grid on
figure(4)
L = length(t_total); 
NFFT = 2^nextpow2(L); 
f = fs*linspace(0,1,NFFT); 

X = fft(OFDM,NFFT)/L; 
S_X = 1./fs*L.*conj((X(1:NFFT)).*conj(X(1:NFFT))); 
S_X = (S_X); 
plot(f/1e9,S_X) 
title([num2str(length(stream)) ' symbols ' num2str(N_channels) ' channels']); 
xlabel('Frequency [GHz]') 
ylabel('S_X(f)') 
grid on


%%
save('myfile.mat','stream','Length','N_channels','x_modulating_with_cpe','cps','x_par_pad','symbol_book','ts')


%% Will be used for the receiver
OFDM_received = OFDM;
% save this information too 
save('OFDM.mat','OFDM')

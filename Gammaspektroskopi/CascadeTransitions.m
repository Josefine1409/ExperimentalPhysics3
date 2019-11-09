clear all; close all; format compact; clc

%% DEMONSTRATE CASCADE TRANSITIONS
%% Import data
[data1, ~, headerlinesOut1] = importdata("manganese/Mn56_???.txt");
[data2, ~, headerlinesOut2] = importdata("manganese/Mn56_???.txt");
time1 = data1.data(:,1);
time2 = data2.data(:,1);
channelCount1 = data1.data(:,2);
channelCount2 = data2.data(:,2);

%% Energy calibration
E1 = @(channel) 0.7364*channel + 0.0625; % keV
uncertantyE1 = ???;
E2 = @(channel) ???*channel + ???; % keV
uncertantyE2 = ???;

%% Energy and time restrictions

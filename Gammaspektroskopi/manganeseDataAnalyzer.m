clear all; close all; format compact; clc

%% Import data and background
[data, ~, headerlinesOut] = importdata("manganese/Mn56_Day1_ch000.txt");
time = data.data(:,1);
channelCount = data.data(:,2);

% background

%% Energy calibration
%

%% Histogram of data with background
figure(1)
hold on
    xlim([0, max(channelCount)])
    histogram(channelCount, 'BinWidth', 1)
hold off

%% Peak fitter
function peakFitter(channelCount, channelToFitMin, channelToFitMax)
    % Only use data between channelToFitMin and channelToFitMax
    channelCountNew = []
    for j = 1:
        if 
            channelCountNew = [channelCountNew; ];
    end
    
    % Fitting function
        % N_f = @(W) alpha * (gamma_f / ((W - M)^2 * c^4 + (gamma^2 / 4)));
        % Background linear
    fun = @(x) c(1) * (c(2) / ((x - c(3))^2 * c^4 + (c(4)^2 / 4)));

    % Inputs to nlinfit
    channelsToFit = channelToFitMin : channelToFitMax;
    guess = ???;

    outputFit = nlinfit(channelsToFit, channelCount, N_f, guess);
end
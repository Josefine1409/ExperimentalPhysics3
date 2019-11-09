clear all; close all; format compact; clc

%% Import data and background
[data, ~, headerlinesOut] = importdata("manganese/Mn56_Day1_ch000.txt");
time = data.data(:,1);
channelCount = data.data(:,2);

%% Energy calibration
E = @(channel) 0.7364*channel + 0.0625; % keV

%% Histogram of data with background
figure(1)
hold on
    xlim([0, max(channelCount)])
    histogram(channelCount, 'BinWidth', 1)
    peakFitter(channelCount, 1130, 1165)
hold off

%% Peak fitter
function peakFitter(channelCount, channelToFitMin, channelToFitMax)
    % Only use data between channelToFitMin and channelToFitMax
    channelCountNew = [];
    for j = 1:length(channelCount)
        if channelCount(j) <= channelToFitMax && channelCount(j) >= channelToFitMin
            channelCountNew = [channelCountNew; channelCount(j)];
        end
    end
    % Sorting the channel counts
    channelCountNewSorted = sort(channelCountNew);
    % Histogram to f(x)
    channelCountNew = accumarray(channelCountNewSorted(:), 1);
    
    % Fitting function
            % N_f = @(W) alpha * (gamma_f / ((W - M)^2 * c^4 + (gamma^2 / 4)));
            % Background linear
        % fun = @(k, C) k(1) * (k(2) / ((C - k(3))^2 + (k(4)^2 / 4))) + k(5)*x + k(6);
        % Linaer + Gauss: a*x + b + a1*exp(-((x-b1)/c1)^2)
    fun = @(k, x) k(4) .* x + k(5) + k(3) .* exp(-((x - k(1)) ./ (k(2))).^2);

    % Inputs to nlinfit
    channelsToFit = channelToFitMin : 1 : channelToFitMax;
    [peakHeight, peakMid] = max(channelCountNew);
    guess = [peakMid, 1, peakHeight, 0, 1]; % [peakMid, peakWidth, peakHeight, slope, offset];

    [values,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(channelsToFit, channelCountNew, fun, guess);
    
    plot(channelsToFit, fun(values, channelsToFit), '-r')
end
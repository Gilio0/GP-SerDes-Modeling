% PAM-4 Signal Generation
clc;
clear;
close all;

% Parameters
numSymbols = 1000; % Number of symbols to generate
symbolRate = 1e3;  % Symbol rate (symbols per second)
Fs = 10e3;         % Sampling frequency (samples per second)

% Generate random data symbols (0, 1, 2, 3)
dataSymbols = randi([0 3], numSymbols, 1);

% Map symbols to PAM-4 levels (-3, -1, 1, 3)
pam4Levels = [-3 -1 1 3];
pam4Signal = pam4Levels(dataSymbols + 1); % Add 1 for index-based mapping

% Upsample the signal for proper sampling rate
samplesPerSymbol = Fs / symbolRate;
pam4SignalUpsampled = upsample(pam4Signal, samplesPerSymbol);

% Apply a rectangular pulse shaping filter
pulseShape = ones(1, samplesPerSymbol); % Rectangular pulse
pam4SignalShaped = conv(pam4SignalUpsampled, pulseShape, 'same');

% Time axis
t = (0:length(pam4SignalShaped)-1) / Fs;

% Plot the PAM-4 Signal in Time Domain
figure;
plot(t, pam4SignalShaped, 'LineWidth', 1.5);
title('PAM-4 Signal in Time Domain');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% Plot Constellation Diagram
figure;
scatter(real(pam4Signal), zeros(size(pam4Signal)), 'filled');
title('PAM-4 Constellation Diagram');
xlabel('In-phase');
ylabel('Quadrature');
grid on;

disp('PAM-4 Signal Generation Complete.')
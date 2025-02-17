%%read the channel and get s21
a='CA_19p75dB_thru.s4p';
b=sparameters(a);
freq=b.Frequencies;
s21 = rfparam(b, 2, 1); % S21: Transmission from Port 1 to Port 2
%% draw the frequence response
figure;
 plot(freq,20*log10(abs(s21)));
 xline(25e9);
 %% get the impulse response
 data_rate = 100e9;
 UI = 2 / data_rate;
    % Define oversample ratio
    samples_per_symbol = 64;
    % Timestep
    dt = UI / samples_per_symbol;
t=0:dt:6000*dt;
 figure;
h=ifft(s21);
h=fftshift(h);
h=real(h);
plot(t,h);
pulse = zeros(1, 6001); % Initialize an array of zeros
pulse2= zeros(1, 6001);
pulse(1:samples_per_symbol) = 3;        % Set the first 64 samples to one
pulse2(samples_per_symbol:samples_per_symbol*2) = 1;        % Set the first 64 samples to one
pulse_response = conv(h, pulse, 'same');
pulse2_response=conv(h, pulse2, 'same');
%% Plot Pulse Response with Cursors
figure;
subplot(2, 1, 1);
plot(t, pulse, 'LineWidth', 1.5);
xline(1e-11);
xlabel('Time (s)');
ylabel('Amplitude');
title('Rectangular Pulse');
grid on;

subplot(2,1,2);
plot(t, pulse_response, 'r', 'LineWidth', 1.5);
hold on;
plot(t, pulse2_response,'b', 'LineWidth', 1.5);
hold off;
xlabel('Time (s)');
ylabel('Amplitude');
title('Pulse Response');
grid on;

%% generate the data
num_bits=1000;
data = randi([0, 1], 1, num_bits); % Generate 10,000 random bits
% Map binary data to 4-PAM levels
scrambler_poly = [1 1 0 0 1]; % Example scrambling polynomial (x^4 + x + 1)
initial_state = [1 0 1 1]; % Initial state of the LFSR

scrambler = comm.Scrambler(2, scrambler_poly, initial_state);
scrambled_bits = scrambler(data.');

% Convert to row vector
data = scrambled_bits.';
data_reshaped = reshape(data, 2, []); % Group bits into pairs
decimal_values = bi2de(data_reshaped', 'left-msb'); % Convert to decimal
pam4_levels = [-3, -1, 1, 3]; % Normalized 4-PAM levels
signal_BR = pam4_levels(decimal_values + 1)'; % Map to PAM levels
signal_BR(1)=3;
signal_BR(2)=-3;
% Data rates in Gbps
data_rate = 100e9;
UI = 2 / data_rate;
    % Define oversample ratio
    fmax=1/2e12;
    samples_per_symbol = 64;

    % Timestep
    dt = UI / samples_per_symbol;
    
    % Oversampled signal
    signal_ideal = repelem(signal_BR, samples_per_symbol); % Oversample signal
    % Plot ideal signal eye diagram
    eyediagram(signal_ideal, samples_per_symbol * 3, UI);
    title(sprintf('%.0fGbps 4-PAM Ideal Signal', data_rate / 1e9));  
    signal_filtered = conv(signal_ideal, h, 'same');
    snr=24;
    signal_noisy = awgn(signal_filtered, snr, 'measured') ;
    signal_noisy=signal_noisy*3/2;
    eyediagram(signal_filtered, samples_per_symbol * 3, UI);
        title(sprintf('%.0fGbps 4-PAM Signal after passing the channel', ...
            data_rate / 1e9));
        
    eyediagram(signal_noisy, samples_per_symbol * 3, UI);
    title(sprintf('%.0fGbps 4-PAM Signal after passing the channel', ...
            data_rate / 1e9));
%%   CTLE implementation
%  fp = 25e9;    % Peak frequency (Hz)
%  fz = 10e9;    % Zero frequency (Hz) - adjust for low-frequency boost
%  zeta = 0.2;   % Damping factor - controls peak sharpness
%  fs = 2 * fp;
%  % Convert to angular frequencies
%  wp = 2 * pi * fp;
%  wz = 2 * pi * fz;
% % 
% % % Define Laplace variable
%  s = tf('s');
% % 
% % % Define CTLE Transfer Function
%  H_ctle = (s + wz) / (s^2 + 2*zeta*wp*s + wp^2);
% % 
% % % Normalize DC Gain
%  k = abs(evalfr(H_ctle, 0)); % Compute DC gain
%  H_ctle = H_ctle / k; % Normalize
% % 
% % % Define frequency vector
%  w = 2 * pi * freq; % Convert to angular frequency
% % 
% % % Compute Frequency Response
%  [H_mag, ~] = freqresp(H_ctle, w); % Get magnitude response
%  H_mag=squeeze(H_mag);
% % get the efffect of it on the channel
%  H_mag_channel=H_mag.*s21;
%  impulse_response_ctle= ifft(H_mag);
%  impulse_response_ctle=fftshift(impulse_response_ctle);
%  impulse_response_ctle=real(impulse_response_ctle);
% % 
% % % Convert to dB
%  H_dB = 20 * log10(squeeze(abs(H_mag))); 
%  impulse_response_after_ctle= ifft(H_mag_channel);
%  impulse_response_after_ctle=fftshift(impulse_response_after_ctle);
%  impulse_response_after_ctle=real(impulse_response_after_ctle);
%  
%  signal_ctle = conv(signal_noisy, impulse_response_ctle, 'same');
%  %signal_ctle=50*signal_ctle;
%  % 
%  
%      eyediagram(signal_ctle, samples_per_symbol * 3, UI);
%   title(sprintf('%.0fGbps 4-PAM Signal after channel and CTLE', ...
%        data_rate / 1e9));
% %plot the frequency response of the CTLE and before it and it's effect
%  figure;
% subplot(2, 1, 1);
%  semilogx(freq, H_dB', 'LineWidth', 1.5);
%  xlabel('Frequency (Hz)');
%  ylabel('Magnitude (dB)');
%  title('CTLE Frequency Response');
%  grid on;
%  xline(25e9, '--r', '25 GHz Peak'); % Mark the peak frequency
% % 
% % % % Plot Impulse Response
% %  subplot(2, 1, 2);
% %  plot(t, impulse_response_ctle, 'r', 'LineWidth', 1.5);
% %  xlabel('Time (s)');
% %  ylabel('Amplitude');
% %  title('CTLE Impulse Response');
% %  grid on;   
% % % 
% figure;
%  subplot(2, 1, 1);
%  plot(freq,20*log10(abs(s21)));
%  xline(25e9);
%   grid on; 
%  subplot(2, 1, 2);
%  plot(freq, 20*log10(abs(H_mag_channel)), 'r', 'LineWidth', 1.5);
%  xline(25e9);
%  title('CTLE effect');
%  grid on; 
 %% sampling the ADC
 %Downsample the signal by selecting the middle sample from each 64-sample block
downsample_factor = 64; % Downsample by a factor of 64
middle_sample_index = downsample_factor / 2; % Middle sample index (32nd sample)

signal_ctle_trimmed = signal_noisy(623:end); % Trim the signal to a multiple of 64

% Extract the middle sample from each 64-sample block
signal_downsampled = signal_ctle_trimmed(middle_sample_index:downsample_factor:end);

% downsample_factor = 64; % Downsample by a factor of 64
% 
% % Trim the signal to a multiple of 64
% signal_ctle_trimmed = signal_noisy(623:end); 
% signal_ctle_trimmed = signal_ctle_trimmed(1:floor(length(signal_ctle_trimmed)/downsample_factor) * downsample_factor);
% 
% % Reshape into 64-sample blocks
% signal_blocks = reshape(signal_ctle_trimmed, downsample_factor, []);
% 
% % Find the index of the max absolute value in each block
% [~, max_idx] = max(abs(signal_blocks), [], 1);
% 
% % Extract the corresponding values from the original blocks
% signal_downsampled = signal_blocks(sub2ind(size(signal_blocks), max_idx, 1:size(signal_blocks,2)));
% 
% % Convert to column vector
% signal_downsampled = signal_downsampled(:);
%% check the symbol error rate

 signal_downsampled_aligned=signal_downsampled;
 signal_BR_aligned=signal_BR;

% Ensure both signals have the same length
min_length = min(length(signal_BR_aligned), length(signal_downsampled_aligned));
signal_BR_aligned = signal_BR_aligned(1:min_length);
signal_downsampled_aligned = signal_downsampled_aligned(1:min_length);
% 
% % Step 2: Quantize the downsampled signal to PAM4 levels
% pam4_levels = [-3, -1, 1, 3]; % PAM4 levels
% signal_downsampled_quantized = zeros(size(signal_downsampled_aligned));
% for i = 1:length(signal_downsampled_aligned)
%     [~, idx] = min(abs(signal_downsampled_aligned(i) - pam4_levels)); % Find the nearest level
%     signal_downsampled_quantized(i) = pam4_levels(idx); % Assign the nearest level
% end
% % Step 3: Compare the quantized downsampled signal to the original signal
% symbol_errors_after_quantization = sum(signal_BR_aligned ~= signal_downsampled_quantized); % Count bit errors
% total_symbols_after_quantization = length(signal_BR_aligned); % Each PAM4 symbol represents 2 bits
% % Step 4: Calculate the Bit Error Rate (BER)
% SER = symbol_errors_after_quantization / total_symbols_after_quantization;
% % Display the results
% disp(['Number of symbols errors: ', num2str(symbol_errors_after_quantization)]);
% disp(['symbol Error Rate (SER): ', num2str(SER)]);
%% ADC after sampling >> quantization
% Step 1: Define ADC parameters
adc_bits = 8; % 8-bit ADC
adc_levels = 2^adc_bits; % 256 levels
adc_min = -3; % Minimum value of the ADC range
adc_max = 3; % Maximum value of the ADC range

% Step 2: Generate ADC levels
adc_step = (adc_max - adc_min) / (adc_levels - 1); % Step size between ADC levels
adc_quantization_levels = linspace(adc_min, adc_max, adc_levels); % ADC levels from -3 to 3

% Step 3: Quantize the downsampled signal to ADC levels
signal_quantized_adc = zeros(size(signal_downsampled_aligned)); % Initialize quantized signal
for i = 1:length(signal_downsampled_aligned)
    % Find the nearest ADC level for each value
    [~, idx] = min(abs(signal_downsampled_aligned(i) - adc_quantization_levels));
    signal_quantized_adc(i) = adc_quantization_levels(idx);
end
   signal_quantized_adc_oversampled = interp(signal_quantized_adc, samples_per_symbol); % Oversample signal

  eyediagram(signal_quantized_adc_oversampled, samples_per_symbol, UI);
  title(sprintf('%.0fGbps 4-PAM Signal after channel and CTLE and ADC', ...
       data_rate / 1e9))

% Step 4: Map the quantized ADC values to the nearest PAM4 level
signal_mapped_pam4 = zeros(size(signal_quantized_adc)); % Initialize mapped PAM4 signal
for i = 1:length(signal_quantized_adc)
    % Find the nearest PAM4 level for each ADC level
    [~, idx] = min(abs(signal_quantized_adc(i) - pam4_levels));
    signal_mapped_pam4(i) = pam4_levels(idx);
end
symbol_errors_after_quantization = sum(signal_BR_aligned ~= signal_mapped_pam4); % Count bit errors
total_symbols_after_quantization = length(signal_BR_aligned); % Each PAM4 symbol represents 2 bits
% Step 4: Calculate the Bit Error Rate (SER)
SER_after_quantization = symbol_errors_after_quantization / total_symbols_after_quantization;
% Display the results
disp(['Number of symbols errors after quantization: ', num2str(symbol_errors_after_quantization)]);
disp(['symbol Error Rate (SER) after quantization: ', num2str(SER_after_quantization)]);
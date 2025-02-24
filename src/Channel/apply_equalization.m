% function signal_ctle = apply_equalization(signal_noisy, config)
%     % CTLE parameters
%     fp = 25e9; % Peak frequency (Hz)
%     fz = 10e9; % Zero frequency (Hz)
%     zeta = 0.2; % Damping factor
% 
%     % CTLE transfer function
%     s = tf('s');
%     H_ctle = (s + 2 * pi * fz) / (s^2 + 2 * zeta * 2 * pi * fp * s + (2 * pi * fp)^2);
%     k = abs(evalfr(H_ctle, 0)); % Normalize DC gain
%     H_ctle = H_ctle / k;
% 
%     % Compute frequency response
%     w = 2 * pi * config.freq;
%     H_mag = squeeze(freqresp(H_ctle, w));
% 
%     % Apply CTLE to channel response
%     H_mag_channel = H_mag .* config.s21;
%     impulse_response_ctle = ifft(H_mag_channel);
%     impulse_response_ctle = fftshift(impulse_response_ctle);
%     impulse_response_ctle = real(impulse_response_ctle);
% 
%     % Apply CTLE to signal
%     signal_ctle = conv(signal_noisy, impulse_response_ctle, 'same');
% end
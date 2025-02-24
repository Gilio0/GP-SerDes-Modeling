function [signal_filtered, signal_noisy,h] = simulate_channel(signal_BR, config)
    % Simulate channel using S-parameters
    h = ifft(config.s21);
    h = fftshift(h);
    h = real(h);

    % Oversample signal
    UI = 2 / config.data_rate;
    dt = UI / config.samples_per_symbol;
    signal_ideal = repelem(signal_BR, config.samples_per_symbol);

    % Apply channel impulse response
    signal_filtered = conv(signal_ideal, h, 'same');

    % Add noise
    signal_noisy = awgn(signal_filtered, config.snr, 'measured');
    signal_noisy = signal_noisy * 3 / 2;
end
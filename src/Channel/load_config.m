function config = load_config()
    % Configuration Structure
    config = struct();
    config.data_rate = 100e9; % Data rate in bps
    config.modulation = 'PAM4'; % Modulation scheme
    config.snr = 20; % SNR in dB
    config.samples_per_symbol = 64; % Oversampling ratio
    config.adc_bits = 8; % ADC resolution
    config.adc_min = -3; % ADC minimum value
    config.adc_max = 3; % ADC maximum value
    config.mu=0.001;
    config.FFE_taps=5;
    config.DFE_taps=2;
    % Load S-parameters
    a = 'CA_19p75dB_thru.s4p';
    b = sparameters(a);
    config.freq = b.Frequencies;
    config.s21 = rfparam(b, 2, 1); % S21 parameter
end
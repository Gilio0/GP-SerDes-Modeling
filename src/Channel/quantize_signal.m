function [signal_quantized_adc, signal_mapped_pam4] = quantize_signal(signal_noisy, config, tb)
    % Downsample the signal
    downsample_factor = config.samples_per_symbol;
    middle_sample_index = downsample_factor / 2;
    signal_noisy_trimmed = signal_noisy(623:end); % Trim the signal
    signal_downsampled = signal_noisy_trimmed(middle_sample_index:downsample_factor:end);

    % Quantize to ADC levels
    adc_levels = 2^config.adc_bits;
    adc_step = (config.adc_max - config.adc_min) / (adc_levels - 1);
    adc_quantization_levels = linspace(config.adc_min, config.adc_max, adc_levels);

    signal_quantized_adc = zeros(size(signal_downsampled));
    for i = 1:length(signal_downsampled)
        [~, idx] = min(abs(signal_downsampled(i) - adc_quantization_levels));
        signal_quantized_adc(i) = adc_quantization_levels(idx);
    end

    % Map to PAM4 levels
    signal_mapped_pam4 = zeros(size(signal_quantized_adc));
    for i = 1:length(signal_quantized_adc)
        [~, idx] = min(abs(signal_quantized_adc(i) - tb.pam4_levels));
        signal_mapped_pam4(i) = tb.pam4_levels(idx);
    end
end
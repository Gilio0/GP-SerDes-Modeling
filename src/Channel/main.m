% Main Script
clear all;
clc;

% Load configuration and testbench
config = load_config();
tb = load_testbench();

% Generate input signal
signal_BR = generate_signal(tb);

% Simulate channel
[signal_filtered, signal_noisy,h] = simulate_channel(signal_BR, config);

% Apply equalization
%signal_ctle = apply_equalization(signal_noisy, config);

% Quantize signal
[signal_quantized_adc, signal_mapped_pam4] = quantize_signal(signal_noisy, config, tb);

% Calculate error rate
[SER , symbol_errors] = calculate_error_rate(signal_BR, signal_mapped_pam4, tb);

% Display results
disp(['Symbol Error Rate (SER) before equalization: ', num2str(SER)]);
disp(['Symbol Errors before equalization : ', num2str(symbol_errors)]);

[des_out2,desicicon_in]=equalization_after_adc(signal_quantized_adc,signal_BR,config,tb);

[SER , symbol_errors] = calculate_error_rate(signal_BR, des_out2, tb);

% Display results
disp(['Symbol Error Rate (SER) after equalization: ', num2str(SER)]);
disp(['Symbol Errors after equalization : ', num2str(symbol_errors)]);


plot_results(repelem(signal_BR, config.samples_per_symbol), signal_filtered, signal_noisy, signal_quantized_adc,desicicon_in, config, tb,h);
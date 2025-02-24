function [SER ,symbol_errors] = calculate_error_rate(signal_BR, signal_mapped_pam4, tb)
    % Align signals
    min_length = min(length(signal_BR), length(signal_mapped_pam4));
    signal_BR_aligned = signal_BR(1:min_length);
    signal_mapped_pam4_aligned = signal_mapped_pam4(1:min_length);

    % Calculate symbol errors
    symbol_errors = sum(signal_BR_aligned ~= signal_mapped_pam4_aligned);
    total_symbols = length(signal_BR_aligned);

    % Calculate SER
    SER = symbol_errors / total_symbols;
end
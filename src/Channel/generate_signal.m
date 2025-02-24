function signal_BR = generate_signal(tb)
    % Generate random bits
    data = randi([0, 1], 1, tb.num_bits);
    
    % Scramble the bits
    scrambler = comm.Scrambler(2, tb.scrambler_poly, tb.initial_state);
    scrambled_bits = scrambler(data.');
    data = scrambled_bits.';
    
    % Map bits to PAM4 levels
    data_reshaped = reshape(data, 2, []);
    decimal_values = bi2de(data_reshaped', 'left-msb');
    signal_BR = tb.pam4_levels(decimal_values + 1)';
end
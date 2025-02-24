function tb = load_testbench()
    % Testbench Structure
    tb = struct();
    tb.num_bits = 50000; % Number of bits to simulate
    tb.num_bits_in_ffe=50000-20;
    tb.seq_length_in_ffe=25000-10;
    tb.scrambler_poly = [1 1 0 0 1]; % Scrambler polynomial
    tb.initial_state = [1 0 1 1]; % Initial state of the scrambler
    tb.pam4_levels = [-3, -1, 1, 3]; % PAM4 levels
end
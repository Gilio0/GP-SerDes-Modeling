% Parameters
% num_symbols = 2; % Number of PAM-4 symbols
% sps = 4; % Samples per symbol
% mu = 0.01; % Step size for phase adjustment
% symbol_rate = 1e6; % Symbol rate (1 Msym/s)
% fs = sps * symbol_rate; % Sampling frequency

% % Generate random PAM-4 symbols
% pam4_symbols = randi([0, 3], 1, num_symbols); % Symbols: 0, 1, 2, 3
% pam4_symbols = 2 * pam4_symbols - 3; % Map to PAM-4 levels: -3, -1, 1, 3

% % Transmit signal (no pulse shaping for wireline)
% tx_signal = kron(pam4_symbols, ones(1,sps )); % Upsample by sps

% % Simulate wireline channel (simple low-pass filter to introduce ISI)
% channel = [1, 0.5, 0.2]; % Simple channel impulse response
% rx_signal = filter(channel, 1, tx_signal);

% % Add noise to the received signal
% SNR = 20; % Signal-to-Noise Ratio in dB
% rx_signal = awgn(rx_signal, SNR, 'measured');

% % Initialize variables for timing recovery
% phase_error = zeros(1, num_symbols);
% sampling_phase = 0; % Initial sampling phase
%sample = zeros(1, num_symbols);


sample = [3 0.4 -1.6];
num_symbols = length(sample);
ref_voltages = [-2.5, -1.5, -0.5, 0.5, 1.5, 2.5]; % Example reference voltages (can be adjusted)
symbols_est = zeros(1, num_symbols); % Estimated symbols
Early = zeros(1, num_symbols);
Late = zeros(1, num_symbols);
errdata = zeros(1, num_symbols);
err_up = zeros(1, num_symbols);
err_low = zeros(1, num_symbols);
A = zeros(1, num_symbols);
B = zeros(1, num_symbols);
C = zeros(1, num_symbols);
mode = zeros(1, num_symbols);

% Mueller-Muller Phase Detector Loop
for n = 1:num_symbols
    % % Interpolate to find the sample at the current sampling phase
    % sample_index = (n-1)*sps + sampling_phase;
    
    % % Ensure the sample index is within valid bounds
    % if sample_index < 1
        % sample_index = 1;
    % elseif sample_index > length(rx_signal)
        % sample_index = length(rx_signal);
    % end
    
    %sample(n) = interp1(1:length(rx_signal), rx_signal, sample_index, 'spline');
    
    % ############ data sampler ############3
    if sample(n) < -2
        symbols_est(n) = -3;
		if sample(n) > ref_voltages(1)
			err_up(n) = 1;
		else 
			err_up(n) = 0;
		end
    elseif sample(n) < 0
        symbols_est(n) = -1;
		if sample(n) > ref_voltages(3)
			err_up(n) = 1;
		else 
			err_up(n) = 0;
		end
		if sample(n) < ref_voltages(2)
			err_low(n) = 1;
		else 
			err_low(n) = 0;
		end
    elseif sample(n) < 2
        symbols_est(n) = 1;
		if sample(n) > ref_voltages(5)
			err_up(n) = 1;
		else 
			err_up(n) = 0;
		end
		if sample(n) < ref_voltages(4)
			err_low(n) = 1;
		else 
			err_low(n) = 0;
		end
    else
        symbols_est(n) = 3;
		if sample(n) < ref_voltages(6)
			err_low(n) = 1;
		else 
			err_low(n) = 0;
		end
    end
	%############### error sampler  ############%
	errdata(n) = xor(err_low(n),err_up(n));
end

for n = 3:num_symbols
%###############waveform selector ############%
	if (sample(n)>sample(n-1))&&(sample(n-1)>sample(n-2))
		mode(n) = 001;
	elseif (sample(n)<sample(n-1))&&(sample(n-1)<sample(n-2))
		mode(n) = 100;
	elseif (sample(n)~=sample(n-1))&&(sample(n-1)==sample(n-2))
		mode(n) = 110;
	elseif (sample(n)==sample(n-1))&&(sample(n-1)~=sample(n-2))
		mode(n) = 011;
	else
		mode(n) = 000;
	end
	
%######## phase detector #############
	switch mode(n)
    case 000 
        A(n) = 0;
        B(n) = 0;
        C(n) = 0;
    case 001 
        A(n) = errdata(n-1);
        B(n) = err_low(n-1);
        C(n) = err_up(n-1);
    case 100 
        A(n) = errdata(n-1);
        B(n) = err_up(n-1);
        C(n) = err_low(n-1);
    case 011 
        if errdata(n-1) == 1
            A(n) = 1;
            B(n) = 0;
            C(n) = 1;
        else
            A(n) = 1;
            B(n) = errdata(n);
            C(n) = errdata(n-1);
        end
    case 110 
        if errdata(n-1) == 1
            A(n) = 1;
            B(n) = 1;
            C(n) = 0;
        else
            A(n) = 1;
            B(n) = errdata(n-1);
            C(n) = errdata(n-2);
        end
    otherwise
        A(n) = 0;
        B(n) = 0;
        C(n) = 0;
end
	Early(n) = A(n)&&B(n);
	Late(n) = A(n)&&C(n);
	
end

figure;
subplot(2,1,1);
plot(Late, 'o');
title('late sampling');
xlabel('Symbol Index');
%ylabel('Phase Error');

subplot(2,1,2);
plot(Early, 'o');
title('Early sampling');
xlabel('Symbol Index');
%ylabel('Phase Error');

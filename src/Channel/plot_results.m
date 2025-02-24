function plot_results(signal_ideal, signal_filtered, signal_noisy, signal_quantized_adc,desicicon_in, config, tb,h)
    % Plot Frequency Response
    figure;
    subplot(2, 1, 1);
    plot(config.freq, 20*log10(abs(config.s21)), 'LineWidth', 1.5);
    xline(25e9, '--r', '25 GHz');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    title('Channel Frequency Response (S21)');
    grid on;

    % Plot Eye Diagrams
    UI = 2 / config.data_rate;
    samples_per_symbol = config.samples_per_symbol;  
 dt = UI / samples_per_symbol;
  t=0:dt:6000*dt;
  pulse = zeros(1, 6001); % Initialize an array of zeros
pulse2= zeros(1, 6001);
pulse(1:samples_per_symbol) = 3;        % Set the first 64 samples to one
pulse2(samples_per_symbol:samples_per_symbol*2) = 1;        % Set the first 64 samples to one
pulse_response = conv(h, pulse, 'same');
pulse2_response=conv(h, pulse2, 'same');
   
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
    

    eyediagram(signal_ideal, samples_per_symbol * 3, UI);
    title('Ideal Signal Eye Diagram');

    eyediagram(signal_filtered, samples_per_symbol * 3, UI);
    title('Filtered Signal Eye Diagram');
         

    eyediagram(signal_noisy, samples_per_symbol * 3, UI);
    title('Noisy Signal Eye Diagram');

    %eyediagram(signal_ctle, samples_per_symbol * 3, UI);
   % title('CTLE Equalized Signal Eye Diagram');
    
   signal_quantized_adc_oversampled = interp(signal_quantized_adc, samples_per_symbol); % Oversample signal
      
    % Plot Quantized Signal Eye Diagram
    eyediagram(signal_quantized_adc_oversampled, samples_per_symbol, UI);
    title('Quantized Signal Eye Diagram');
    
    desicicon_in2 = interp(desicicon_in, samples_per_symbol); % Oversample signal

  eyediagram(desicicon_in2, samples_per_symbol, UI);
  title('equalized Signal Eye Diagram');
end
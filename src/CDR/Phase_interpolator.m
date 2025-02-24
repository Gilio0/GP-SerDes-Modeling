clear;
clc;
% Parameters
f = 1; % Frequency of the square wave (Hz)
t = 0:0.001:2; % Time vector from 0 to 2 seconds with 1 ms steps
Early = [1 1 1 1 1 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 0 0 0 0 0];
Late =  [0 0 0 0 0 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
num_symbols = length(Early);
counter = 0;
A = zeros(1, num_symbols);
B = zeros(1, num_symbols);

% Phase interpolator controller 
for i = 1:num_symbols
    if Early(i) == 1
        if counter >= 127
            counter = 0;
        else
            counter = counter + 1;
        end
    elseif Late(i) == 1
        if counter <= 0
            counter = 127;
        else
            counter = counter - 1;
        end
    else 
        counter = counter;
    end
    
    % Phase interpolator 128 steps resolution 
    if (0 <= counter) && (counter <= 8)
        A(i) = 3 * counter;
        B(i) = 64 - 2 * counter;
    elseif (9 <= counter) && (counter <= 16)
        A(i) = 16 + counter;
        B(i) = 64 - 2 * counter;
    elseif (17 <= counter) && (counter <= 24)
        A(i) = 2 * counter;
        B(i) = 48 - counter;
    elseif (25 <= counter) && (counter <= 32)
        A(i) = 2 * counter;
        B(i) = 96 - 3 * counter;
    elseif (33 <= counter) && (counter <= 40)
        A(i) = 128 - 2 * counter;
        B(i) = 96 - 3 * counter;
    elseif (41 <= counter) && (counter <= 48)
        A(i) = 128 - 2 * counter;
        B(i) = 16 - counter;
    elseif (49 <= counter) && (counter <= 56)
        A(i) = 80 - counter;
        B(i) = 64 - 2 * counter;
    elseif (57 <= counter) && (counter <= 64)
        A(i) = 192 - 3 * counter;
        B(i) = 64 - 2 * counter;
    elseif (65 <= counter) && (counter <= 72)
        A(i) = 192 - 3 * counter;
        B(i) = 2 * counter - 192;
    elseif (73 <= counter) && (counter <= 80)
        A(i) = 48 - counter;
        B(i) = 2 * counter - 192;
    elseif (81 <= counter) && (counter <= 88)
        A(i) = 128 - 2 * counter;
        B(i) = counter - 112;
    elseif (89 <= counter) && (counter <= 96)
        A(i) = 128 - 2 * counter;
        B(i) = 3 * counter - 288;
    elseif (97 <= counter) && (counter <= 104)
        A(i) = 2 * counter - 256;
        B(i) = 3 * counter - 288;
    elseif (105 <= counter) && (counter <= 112)
        A(i) = 2 * counter - 256;
        B(i) = counter - 80;
    elseif (113 <= counter) && (counter <= 120)
        A(i) = counter - 144;
        B(i) = 2 * counter - 192;
    elseif (121 <= counter) && (counter <= 128)
        A(i) = 3 * counter - 384;
        B(i) = 2 * counter - 192;
    else 
        A(i) = 0;
        B(i) = 0;
    end

    % Generate the square wave with phase shift
    square_wave = square(2 * pi * f * t - atan2(A(i), B(i)));

    % Plot the square wave for the current iteration
    plot(t, square_wave);
    ylim([-1.5 1.5]); % Adjust y-axis limits for better visualization
    xlabel('Time (s)');
    ylabel('Amplitude');
    title(['Phase Interpolator Output - Iteration: ', num2str(i)]);
    grid on;
    pause(1); % Pause for 1 second to observe the effect
end
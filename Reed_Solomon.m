clc
clear all;
%% main
fcr=0;
 generator=2;
 n=255;
 k=223;
 prim_poly=285;
 m=8;
 %msg_in = randi([0 255], 1, k);
 msg_in =[163   227   195   124    17   132   103   148   247   246   255   131   184    37    47    60   235 ...
     20   178    48   209    83   204   208    81    52    41   233   181    30   163    38    20   190 ...
     49     0     5   102   151    12     0   147    89    41   159   141   217   102    13   255    11 ...
     32   117   205   103   193    86    88   156   240    45    30   117   249   203    40   202   183 ...
      78    36   254   197   131    96    68   152    27   174    52   227   158     9   217    35    59 ...
      96    51   170    32   214   179    75   192   213    69   134     2   190   247   249   169   107 ...
      2    92   176   103    60   244    63     9    67   206   203    32    32    16   212   150    86 ...
      193   173    69   157   121    51   150   220   189   217   131    87    46    20   118    23   150 ...
      215    10    98   135   165    27   193   100     9    16     2   147   194    22    16   163   231 ...
      104   200   129    67    97   120   184   148   119    66   164     9    15   189     5   166   227 ...
      230    75   208   145   238   204    70   255   186    73   151   191   115    10   244   123    21 ...
      13   161    78   189   144   182   165    78   139   204   173   249   169    63   115    64   194 ...
      133    27   175   169   231   125   118   169   240    83     5   197    21   137   204   140   205 ...
      80   101];
 nsym = n-k; % Number of error correction symbols
 max_errors=floor(nsym/2);
% % Encode message using Reed-Solomon encoding
 
encoded_msg = rs_encode_msg(msg_in, nsym,prim_poly, fcr, generator);
%Toolbox Encoder
rsEncoder=comm.RSEncoder('CodewordLength', n, 'MessageLength', k, ...
    'BitInput', false, 'PrimitivePolynomial', [1 0 0 0 1 1 1 0 1], 'PrimitivePolynomialSource', 'Property');
encoded=step(rsEncoder, msg_in.');
if (encoded_msg==encoded.')
    disp("encoding is successful");
else
        disp("encoding is failed");
end
 corrupted=encoded_msg;
 corrupted(10)=15;
%  S = rs_calc_syndromes(corrupted, 30, 0, 2);
%     S=trim(S);
%  
%      g = euclidean_algorithm(S, max_errors);

[decoded, error_pos, error_mag, g, S] = rsDecoder(corrupted, n, k);
if (decoded==msg_in)
    disp("Decoding is successful");
else
        disp("Decoding is failed");
end
%%
function result = gf_add(x, y)
    % GF addition is XOR
    result = bitxor(x, y);
end

function result = gf_multiply(a, b, prim_poly)
    % Multiply two numbers in GF(256) modulo a primitive polynomial
    % a and b are the input numbers (polynomials represented as integers)
    % prim_poly is the primitive polynomial (default is x^8 + x^4 + x^3 + x + 1, represented as 0x11D)

    % Default primitive polynomial for GF(256)
    if nargin < 3
        prim_poly = 285; % 0x11D in decimal
    end

    % Initialize the result
    result = 0;

    % Perform carryless multiplication
    for i = 1:8
        if bitand(b, 1) % If the least significant bit of b is set
            result = bitxor(result, a); % Add a to the result
        end
        a = bitshift(a, 1); % Multiply a by x
        if a >= 256 % If a overflows 8 bits
            a = bitxor(a, prim_poly); % Reduce modulo the primitive polynomial
        end
        b = bitshift(b, -1); % Divide b by x
    end
end

function r = gf_poly_add(p, q)
    % Galois Field polynomial addition (mod 2)
    
    % Ensure r is the correct size, maximum length of p and q
    r = zeros(1, max(length(p), length(q)), 'uint8');
    
    % Copy p into the end of r
    r(end-length(p)+1:end) = p;
    
    % XOR elements of q with r
    for i = 1:length(q)
        r(i + length(r) - length(q)) = bitxor(r(i + length(r) - length(q)), q(i));
    end
end

function r = gf_poly_mul(p, q)
    % Multiply two polynomials in Galois Field
   
    % Pre-allocate the result array
    r = zeros(1, length(p) + length(q)-1);
    
    % Precompute the logarithm of p
    %lp = gf_log(p); % Directly index the lookup table
    
    % Compute the polynomial multiplication
    for j = 1:length(q)
        qj = q(j); % Load the coefficient once
        if qj ~= 0 % log(0) is undefined, check that
            %lq = gf_log(qj); % Precompute log of qj
            for i = 1:length(p)
                if p(i) ~= 0 % log(0) is undefined, check that
                    r(i + j-1) = bitxor(r(i + j-1), gf_multiply(p(i),qj,285));
                end
            end
        end
    end
    r = trim(r);
end

function [quotient, remainder] = gf_poly_div(dividend, divisor)
    % Fast polynomial division using Extended Synthetic Division optimized for GF(2^p)
% Default primitive polynomial for GF(256)
    if nargin < 3
        prim_poly = 285; % 0x11D in decimal
    end
    % Copy the dividend and pad with zeros where the ECC bytes will be computed
    msg_out = dividend;
    
    for i = 1:length(dividend) - (length(divisor) - 1)
        coef = msg_out(i); % Pre-caching
        if coef ~= 0 % log(0) is undefined, avoid it explicitly
            for j = 2:length(divisor) % Skip first coefficient of divisor
                if divisor(j) ~= 0 % log(0) is undefined
                    msg_out(i + j - 1) = bitxor(msg_out(i + j - 1), gf_multiply(divisor(j), coef,prim_poly));
                end
            end
        end
    end
    % Compute the index where quotient and remainder are separated
    separator = length(dividend) - (length(divisor) - 1);
    quotient = trim(msg_out(1:separator));
    remainder = msg_out(separator + 1:end);
end
function result = gf_pow(a, k, prim_poly)
    % Compute a^k in GF(256) modulo a primitive polynomial
    % a is the input polynomial (represented as an integer)
    % k is the exponent
    % prim_poly is the primitive polynomial (default is x^8 + x^4 + x^3 + x + 1, represented as 0x11D)

    % Default primitive polynomial for GF(256)
    if nargin < 3
        prim_poly = 285; % 0x11D in decimal
    end

    % Handle the case where k = 0
    if k == 0
        result = 1; % Any non-zero element to the power of 0 is 1 in GF(256)
        return;
    end

    % Initialize the result
    result = 1;

    % Perform exponentiation by squaring
    while k > 0
        % If the current bit of k is set, multiply the result by a
        if bitand(k, 1)
            result = gf_multiply(result, a, prim_poly);
        end

        % Square a
        a = gf_multiply(a, a, prim_poly);

        % Shift k to the right
        k = bitshift(k, -1);
    end
end
function y = gf_poly_eval(poly, x)
    % Evaluates a polynomial in GF(2^p) using Horner's scheme    
        if nargin < 3
        prim_poly = 285; % 0x11D in decimal
    end

    y = poly(1);
    for i = 2:length(poly)
        y = bitxor(gf_multiply(y, x,prim_poly), poly(i));
    end
end

%% encoding
function g = rs_generator_poly(nsym, fcr, generator)
    % Generate an irreducible generator polynomial for Reed-Solomon encoding

    if nargin < 2
        fcr = 0;
    end
    if nargin < 3
        generator = 2;
    end
    
    g = [1];
    for i = 1:(nsym)
        g = gf_poly_mul(g, [1, gf_pow(generator, i + fcr,285)]);
    end
end

function msg_out = rs_encode_msg(msg_in, nsym,prim_poly, fcr, generator)
    % Reed-Solomon encoding using polynomial division (Extended Synthetic Division)  

    if nargin < 4
        fcr = 0;
    end
    if nargin < 5
        generator = 2;
    end
    gen = rs_generator_poly(nsym, fcr, generator);
        
    msg_out = [msg_in, zeros(1, length(gen)-1)]; % Initialize msg_out with msg_in and pad with len(gen)-1 zeros
    disp(gen);
    % Extended synthetic division main loop
    for i = 1:length(msg_in)
        coef = msg_out(i); % Use msg_out for updated values
        if coef ~= 0
            for j = 2:length(gen) % Skip first coefficient of divisor
     msg_out(i + j - 1) = bitxor(msg_out(i + j - 1), gf_multiply(coef , gen(j),prim_poly));
            end
        end
    end
    % Overwrite the quotient part with the original message bytes
    msg_out(1:length(msg_in)) = msg_in;
end

%% decoding
function syndromes = rs_calc_syndromes(msg, nsym, fcr, generator)
syndromes = zeros(1, nsym);
    for i = 1:nsym
        syndromes(i) = gf_poly_eval(msg, gf_pow(2, i ,285));
    end
    syndromes=trim(syndromes);
end

function result = trim(poly)
    % Remove leading zeros from a polynomial
    first_nonzero = find(poly ~= 0, 1); 
    if isempty(first_nonzero)
    result = [0];  % Ensuring it's a row vector
else
    result = poly(first_nonzero:end);
end

end

% Add leading zeros
function result = pad(poly, target_length)
    % Ensure the polynomial has at least target_length coefficients
    current_length = length(poly);
    if current_length < target_length
        result = [zeros(1, target_length - current_length), poly];
    else
        result = poly;
    end
end

function g = euclidean_algorithm(S, max_errors)
    % Euclidean algorithm to compute the error locator polynomial g(x)
    % S: Syndrome polynomial (vector of syndromes)
    % max_errors: Maximum number of errors the code can correct

    % Ensure S is not empty
    if isempty(S) || all(S == 0)
        disp('No errors detected, returning g(x) = 1.');
        g = [1]; % No errors, return g(x) = 1
        return;
    end

    % Initialize polynomials
    r0 = [1, zeros(1, 2 * max_errors)]; % r0(x) = x^(2t)
    r0 = trim(r0);
    r1 = trim(S); % Syndrome polynomial
    size_r0 = length(r0);
    
    g0 = zeros(1, size_r0); % g0(x) = 0
    g1 = [zeros(1, size_r0 - 1), 1]; % g1(x) = 1

    % Iterate until stopping condition is met
    while true
        % Ensure r1 is valid for division
        if r1(1) == 0
            error('Leading coefficient of r1 is zero, division undefined!');
        end
        
        % Perform polynomial long division in GF(256)
        [quotient, remainder] = gf_poly_div(r0, r1);
        quotient = pad(quotient, length(g1));
        
        % Compute new error locator polynomial
        c = gf_poly_mul(quotient, g1);
        c = trim(c);
        c = pad(c, length(g0));
        g = gf_poly_add(g0, c);
disp(['Quotient: ', num2str(quotient)]);
disp(['g(x) before update: ', num2str(g0)]);
disp(['g(x) after update: ', num2str(g)]);

        % Stop when degree of remainder is < max_errors
        if length(remainder) < max_errors || all(remainder == 0)
    break;
end



        % Update values for next iteration
        r0 = trim(r1);
        r1 = trim(remainder);
        g0 = g1;
        g1 = g;
    end

    % Trim leading zeros from the result
    g = trim(g);
if length(g) > max_errors
    error('Error locator polynomial degree exceeds max_errors.');
end
test_dividend = [1 2 3 4 5];
test_divisor = [1 1];
[quotient, remainder] = gf_poly_div(test_dividend, test_divisor);
disp(['Quotient: ', num2str(quotient)]);
disp(['Remainder: ', num2str(remainder)]);

disp(['Current remainder: ', num2str(remainder)]);
disp(['Current g(x): ', num2str(g)]);

    % Ensure g(x) is valid
    if isempty(g) || length(g) < 2
        error('Euclidean algorithm output polynomial g(x) is too short!');
    end
end

%% Find Zeros of Polynomial
 function error_pos = find_zeros(g)
    error_pos = [];
for i = 255:-1:1
    if gf_poly_eval(g, gf_pow(2, i)) == 0
        error_pos = [error_pos, i]; % Append instead of overwriting
    end
end

    disp(error_pos);
end

%% Solve Error Magnitudes
 function error_mag = solve_error_magnitudes(S, error_pos)
    size_error = length(error_pos);
    syndrome_vals = S;
    b(:, 1) = syndrome_vals(1:size_error);
    for i = 1:size_error
    E = 2.^(i * (255 - error_pos));
    E(i, :) = E;
end
error_mag = E \ b';

 end
function [decoded, error_pos, error_mag, g, S] = rsDecoder(encoded, n, k)
    
    % Initialize Galois Field
    
    max_errors = floor((n - k) / 2);
    errors = zeros(1, n);
    g=[];
    % Compute syndromes
    S = rs_calc_syndromes(encoded, 32, 0, 2);
    S=trim(S);
    if all(S == 0)
        decoded = encoded(1:k);
        error_pos = [];
        error_mag = [];
        return;
    end
    disp('Syndrome Polynomial S:'), disp(S);
    S=trim(S);
    % Compute error locator polynomial using Berlekamp-Massey
    g = euclidean_algorithm(S, max_errors);
    disp('Error Locator Polynomial g:'), disp(g);
    
    % Find error positions
    error_pos = find_zeros(g);
    if isempty(error_pos)
        decoded = encoded(1:k);
        error_mag = [];
        return;
    end
    
    % Solve for error magnitudes
    error_mag = solve_error_magnitudes(S, error_pos);
    
    % Correct errors
    display(error_pos);
    errors(error_pos) = error_mag;
    errors_full = zeros(1, length(encoded));  % Initialize a zero vector
    errors_full(error_pos) = error_mag;       % Assign magnitudes to error positions
    corrected = bitxor(encoded, errors_full); % Apply correction

    decoded = corrected(1:k);
end

function d = poly_degree(p)
    % Find the highest nonzero term index in the polynomial p
    idx = find(p ~= 0, 1, 'last'); % Get the last nonzero coefficient index
    if isempty(idx)
        d = -1; % Define degree of zero polynomial as -1
    else
        d = idx - 1; % MATLAB indices start from 1, so subtract 1
    end
end


%generation data stream
num_bits =50000-20 ;
binary_data_stream = randi( [0,1] , 1,num_bits );
binary_data_stream2 = randi( [0,1] , 1,num_bits );
%% parameter
Eb=1;
%snr_db= -5 : 1 : 20 ;
snr_db=5;         %no noise
snr_linear = 10 .^( snr_db / 10 );
N0 = Eb ./snr_linear;
mu=0.001;
%% pam4 mapping
modulated_data=PAM4_mapping(binary_data_stream);
modulated_data2=PAM4_mapping(binary_data_stream2);
pam4_table = [(-3), (-1), (1) , (3)];
%% filter taps
FFE_taps=5;
DFE_taps=2;
%% channel response
h_channel=[1.2 0.7 -0.5];
channel_length=length(h_channel);
channel_out=filter(h_channel,1,modulated_data);
channel_out2=filter(h_channel,1,modulated_data2);
%%adding AWGN noise
Es=( sum( ( real(modulated_data) .^ 2 ) )  ) / (num_bits/2) ;
E_avg=Es/2;
noise_inphase = randn(length(modulated_data),1) .* sqrt( (N0/2)*E_avg  )   ;
channel_out_noise=channel_out+noise_inphase;

Es=( sum( ( real(modulated_data2) .^ 2 ) )  ) / (num_bits/2) ;
E_avg=Es/2;
noise_inphase = randn(length(modulated_data2),1) .* sqrt( (N0/2)*E_avg  )   ;
channel_out_noise2=channel_out2+noise_inphase;

%% channel output
num_data=length(modulated_data);
FFE_output=zeros(1,num_data);
DFE_output=zeros(1,num_data);
desicicon_in=zeros(1,num_data);
des_out = zeros(1,num_data);
error = zeros(1,num_data);
% ffe_history=zeros(num_data,FFE_taps);
% dfe_history=zeros(num_data,DFE_taps);

%% training seq
seq_length=24990;
ffe_window=zeros(1,FFE_taps);dfe_window=zeros(1,DFE_taps);
ffe_taps_initial=zeros(FFE_taps,num_data+1);ffe_taps_initial(3,:)=1;
dfe_taps_initial=zeros(DFE_taps,num_data+1);
for n=1:1:seq_length
    ffe_window = [signal_quantized_adc(n), ffe_window(1:end-1)];
    FFE_output(n)= ffe_taps_initial(:,n)' * ffe_window' ;
    DFE_output(n)= dfe_taps_initial(:,n)'* dfe_window' ;
   
    desicicon_in(n)=FFE_output(n)+DFE_output(n);
    
    distance = abs( desicicon_in( n ) - pam4_table );
    [min_value , index ] = min(distance);
    des_out(n)= pam4_table(index);
    
    error(n)= signal_BR(n)- desicicon_in(n);
    
    ffe_taps_initial(:,n+1) = ffe_taps_initial(:,n) + (mu*error(n)) .* ffe_window';
    dfe_taps_initial(:,n+1) = dfe_taps_initial(:,n) + (mu*error(n)) .* dfe_window';
    
    dfe_window = [signal_BR(n), dfe_window(1:end-1)];
end
%%exact
%w_mmse=MMSE_exact(channel_out_noise,modulated_data,21);


%% data
ffe_taps=ffe_taps_initial(:,seq_length);
dfe_taps=dfe_taps_initial(:,seq_length);
des_out2 = zeros(1,num_data);
% ffe_history=zeros(FFE_taps,num_data+1);
% dfe_history=zeros(DFE_taps,num_data+1);
for n=1:1:num_data
    ffe_window = [signal_quantized_adc(n), ffe_window(1:end-1)];
    FFE_output(n)= ffe_taps' * ffe_window' ;
    DFE_output(n)= dfe_taps'* dfe_window' ;
   
    desicicon_in(n)=FFE_output(n)+DFE_output(n);
    
    distance = abs( desicicon_in( n ) - pam4_table );
    [min_value , index ] = min(distance);
    des_out2(n)= pam4_table(index);
    
    error(n)= des_out2(n)- desicicon_in(n);
    
    ffe_taps = ffe_taps + (mu*error(n)) .* ffe_window';
    dfe_taps = dfe_taps + (mu*error(n)) .* dfe_window';
    
    dfe_window = [des_out2(n), dfe_window(1:end-1)];
end
%%
num_error=sum(des_out2 ~=signal_BR(1:end-10)' );
%disp(num_error);
fprintf('number of error = %d\n',num_error);
%%
figure('name','Equalizers taps' );
x_axis=1:num_data+1;
subplot(2,1,1);
plot(x_axis,ffe_taps_initial);
title('FFE taps');
ylabel('taps weight' );
xlabel('iterations');
grid on ;

subplot(2,1,2);
plot(x_axis,dfe_taps_initial);
title('DFE taps');
ylabel('taps weight' );
xlabel('iterations');
grid on ;

   desicicon_in2 = interp(desicicon_in, samples_per_symbol); % Oversample signal

  eyediagram(desicicon_in2, samples_per_symbol, UI);
  title(sprintf('%.0fGbps 4-PAM Signal after channel and CTLE and ADC and DSP', ...
       data_rate / 1e9))
   
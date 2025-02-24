function [des_out2,desicicon_in]=equalization_after_adc(signal_quantized_adc,signal_BR,config,tb)
mu=config.mu;
FFE_taps=config.FFE_taps;
DFE_taps=config.DFE_taps;
num_data=tb.seq_length_in_ffe;
seq_length=tb.seq_length_in_ffe;
pam4_table= tb.pam4_levels;
FFE_output=zeros(1,num_data);
DFE_output=zeros(1,num_data);
desicicon_in=zeros(1,num_data);
des_out = zeros(1,num_data);
error = zeros(1,num_data);

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

%% data
ffe_taps=ffe_taps_initial(:,seq_length);
dfe_taps=dfe_taps_initial(:,seq_length);
des_out2 = zeros(1,num_data);

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
des_out2=des_out2';
end
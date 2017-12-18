clc
clear
close all
% begin
% N_packet = 1000; % No of iterations
N_packet = 1000;
b = 4; % modulation index 1:BPSK, 2:QPSK, 4: 16 QAM, 6: 64 QAM
N_frame = 4; % No of Modulation symbols per packet
M = 16;
SNRdBs = (0:5:20);
sq2 = sqrt(2);
mode = 'sdm';
channel = 'rayleigh';
plotconst = 'on';

for i_SNR = 1:length(SNRdBs)
    SNRdB = SNRdBs(i_SNR);
    sigma = sqrt(0.5/(10^(SNRdB/10)));
    for i_packet = 1:N_packet
       % Transmitter
       % build msg_symbol as bits block
       % 1     0     0     1
       % 0     0     1     0
       
       msg_symbol = randi([0 1],N_frame*b,1);
       conv_encoded_bits = convencode(msg_symbol');
       
       tx_bits = reshape(conv_encoded_bits,[],4);
       tx_bits = tx_bits';
       %%%%%%%%%%%%%%%%%%%
       % NO ENVIRONMENT
       %%%%%%%%%%%%%%%%%%%
       
       % re-organization bit structure change to properly format to send to
       % viterbi function.
       new_soft_rx_bits =[];
       for i=1:1:size(tx_bits)
           new_soft_rx_bits = [new_soft_rx_bits tx_bits(i,:)];
       end
       
       % change from 1     1     1     1     1     1     0     0 to
       % 1     1     1     0
       % 1     1     1     0
       odd_new_soft_rx_bits = new_soft_rx_bits(:,1:2:end);  % odd matrix
       even_new_soft_rx_bits = new_soft_rx_bits(:,2:2:end);  % event matrix
       new_soft_rx_bits_2 = [odd_new_soft_rx_bits ; even_new_soft_rx_bits];
        
       viterbi_bits = myviterbi(new_soft_rx_bits_2');
       %keep error of each i packet, then we calculate BER
       soft_error = sum(abs(msg_symbol' - viterbi_bits))/2;
       
       %hard_errors(i_packet) = sum(hard_error);
       soft_errors(i_packet) = sum(soft_error);
    end %end for loop for i packet
    
    soft_BER(i_SNR) = sum(soft_errors)/(N_packet*N_frame*b*2*2);
end %end for loop for i SNR

figure
hold on
%Print SOFT decision BER vs SNR
semilogy(SNRdBs, soft_BER, 'r^-')
grid on;
xlabel('SNR[dB]');
ylabel('BER');
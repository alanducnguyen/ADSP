clc
clear
close all
% begin
% N_packet = 1000; % No of iterations
N_packet = 1000;
b = 4; % modulation index 1:BPSK, 2:QPSK, 4: 16 QAM, 6: 64 QAM
N_frame = 2; % No of Modulation symbols per packet
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
       tx_sym = mapper(tx_bits');
       
       X =  zeros(N_frame*b*2,1);
       % Transmiter signal
       X = tx_sym;
       
       %channel
       H = zeros(N_frame*2,1); % RX antenal channel
       
       if strcmp(channel,'rayleigh')
           H(:,1) = (randn(N_frame*2,1)+ 1i*randn(N_frame*2,1))/sq2;
       elseif strcmp(channel,'awgn')
           H = repmat([1,0],N_frame,1);
       else
           error('This channel is not supported');
       end
       
       %Receiver
       R =  H.*X+sigma*(randn(length(X),1)+1i*randn(length(X),1)); %noise = signma*awgn

       %Retrieve rx signal with noise
       rx_sym = R./H;
       
       %hard_rx_bits = harddemapper(rx_sym);
       soft_rx_bits = softdemapper(rx_sym);
       
       % re-organization bit structure change to properly format to send to
       % viterbi function.
       new_soft_rx_bits =[];
       for i=1:1:size(soft_rx_bits)
           new_soft_rx_bits = [new_soft_rx_bits soft_rx_bits(i,:)];
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
clc
clear
close all
% begin
% N_packet = 1000; % No of iterations
N_packet = 1000;
b = 4; % modulation index 1:BPSK, 2:QPSK, 4: 16 QAM, 6: 64 QAM
N_frame = 2; % No of Modulation symbols per packet
M = 16;
SNRdBs = (10:5:40);
sq2 = sqrt(2);
mode = 'sdm';
channel = 'rayleigh';
plotconst = 'on';
is_convolution_encoding = "on";

for i_SNR = 1:length(SNRdBs)
    SNRdB = SNRdBs(i_SNR);
    sigma = sqrt(0.5/(10^(SNRdB/10)));
    for i_packet = 1:N_packet
       % Transmitter
       % build msg_symbol as bits block
       % 1     0     0     1
       % 0     0     1     0
       
       %msg_symbol = randi([0 1],N_frame*b,1);
%        msg_symbol = round(rand([1 N_frame*b]));
%        
%        conv_encoded_bits = convencode(msg_symbol);
%        tx_bits = reshape(conv_encoded_bits,[4 64]);
%        
%        tx_sym = mapper(tx_bits');

       tx_error_symbols = [];
       hard_error_symbols = [];
       
       msg_symbol = round(rand([1 N_frame*b]));
       
       conv_encoded_bits = convencode(msg_symbol);
       
       % conv_encoded_bits = [0
       %                      1
       %                      0
       %                      1]
       %
       tx_sym = mapper(conv_encoded_bits,b,N_frame,is_convolution_encoding);
       
       X =  zeros(N_frame*b,1);
       % Transmiter signal
       X = tx_sym;
       
       %channel
       H = zeros(N_frame*2,1); % RX antenal channel
       
       if strcmp(channel,'rayleigh')
           H(:,1) = (randn(N_frame*2,1)+ 1i*randn(N_frame*2,1))/sq2;
       elseif strcmp(channel,'awgn')
           H = repmat([1,0],N_frame*2,1);
       else
           error('This channel is not supported');
       end
       
       %Receiver
       N =  sigma*(randn(length(X),1)+randn(length(X),1)*1i); %noise = signma*awgn
       
       R =  H.*X+N; 

       %Retrieve rx signal with noise
       rx_sym = R./H;
       
       %hard_rx_bits = harddemapper(rx_sym);
       soft_rx_bits = softdemapper(rx_sym);
       soft_error_mapper = sum(abs(conv_encoded_bits - soft_rx_bits));
       soft_error_mappers(i_packet) = sum(soft_error_mapper);
       % re-organization bit structure change to properly format to send to
       % viterbi function.
       soft_rx_bits_reshape = reshape(soft_rx_bits,[2,length(soft_rx_bits)/2]);

        % viterbi_bits =
        % 0     1     0     1     1     1     0     1     0     1     0     0     1     1     0     0
       viterbi_bits = myviterbi(soft_rx_bits_reshape.');

       %keep error of each i packet, then we calculate BER
       soft_error_conv = sum(abs(msg_symbol - viterbi_bits));
       if soft_error_conv>0
           msg_symbol;
           viterbi_bits;
       end
       soft_errors_conv(i_packet) = sum(soft_error_conv);
    end %end for loop for i packet
    soft_mapper_BER(i_SNR) = sum(soft_error_mappers)/(N_packet*N_frame*b*2);
    soft_BER_conv(i_SNR) = sum(soft_errors_conv)/(N_packet*N_frame*b);
    
end %end for loop for i SNR


axes2 = axes('Parent',figure,'Yscale','log','YMinorTick','on','YminorGrid','on','FontSize',12,'FontName','Time New Roman');
xlim(axes2,[min(SNRdBs) max(SNRdBs)]);
ylim(axes2,[1e-006 1]);
grid(axes2,'on');
hold(axes2,'on');
%Print BER after 16 QAM soft demapper
semilogy(SNRdBs, soft_mapper_BER, 'bo-')
hold on;
%Print BER after convolution-viterbi
semilogy(SNRdBs, soft_BER_conv, 'ro-')
hold on;
grid on;
xlabel('SNR[dB]');
ylabel('BER');
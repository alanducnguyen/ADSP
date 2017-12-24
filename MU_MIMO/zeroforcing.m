clc
clear
close all
% begin
% N_packet = 1000; % No of iterations
N_packet = 3000;
b = 4; % modulation index 1:BPSK, 2:QPSK, 4: 16 QAM, 6: 64 QAM
N_frame = 1; % No of Modulation symbols per packet
M = 16;
SNRdBs = (0:5:40);
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

       tx_error_symbols = [];
       hard_error_symbols = [];
       
       msg_symbol = round(rand([2 N_frame*b]));
       
       [row, col] = size(msg_symbol);
       mimo_conv_encoded_bits=zeros(col*2,row);
       
       for i=1:1:row
           conv_encoded_bits = convencode(msg_symbol(i,:));
           mimo_conv_encoded_bits(:,i)= conv_encoded_bits;
       end
       % conv_encoded_bits = [0
       %                      1
       %                      0
       %                      1]
       mimo_tx_syms = zeros(N_frame*2,2);
       for idx=1:row
            % tx_sym =
            % -3.0000 + 1.0000i
            % 3.0000 + 1.0000i
            % 1.0000 - 1.0000i
            % -1.0000 - 3.0000i
            % -1.0000 + 3.0000i
            % -1.0000 - 3.0000i
            % 3.0000 + 1.0000i
            % -3.0000 + 1.0000i
            
           tx_sym = mapper(mimo_conv_encoded_bits(:,idx),b,N_frame);
           mimo_tx_syms(:,idx) = tx_sym;
       end
       
        % mimo_tx_syms =
        %       X1                   X2
        % -3.0000 - 3.0000i  -3.0000 - 3.0000i
        % 1.0000 - 1.0000i   1.0000 + 3.0000i
       % Transmiter signal
       % X11 X12
       % X21 X22
       mimo_tx_syms = mimo_tx_syms.';
       
       X = mimo_tx_syms;
       
       %channel
       H = zeros(N_frame*2,2); % RX antenal channel
       
       if strcmp(channel,'rayleigh')
           H = (randn(N_frame*2,2)+ 1i*randn(N_frame*2,2))/sq2;
       elseif strcmp(channel,'awgn')
           H = repmat([1,0],N_frame*2,2);
       else
           error('This channel is not supported');
       end
       %noise = signma*awgn
       N = sigma*(randn(2,1)+randn(2,1)*1i);
       %Receiver
       R =  H*X+N; 
       
       H_inv = H^-1;
       mimo_rx_sym = H_inv*R;
       
       %hard_rx_bits = harddemapper(rx_sym);
        % rx_sym =
        % -2.1500 - 3.6662i
        % -3.6876 + 0.4974i
        % 2.9912 - 5.3461i
        % -0.3573 + 3.2709i
        % -3.3309 + 2.1055i
        % 0.3879 - 1.0064i
        % 1.4661 - 1.5155i
        % -2.0477 + 4.0303i
       mimo_soft_rx_bits = zeros(N_frame*b*2,2);
       for rx_idx=1:1:row
           rx_sym = mimo_rx_sym(rx_idx,:);
           soft_rx_bits = softdemapper(rx_sym.');
           mimo_soft_rx_bits(:,rx_idx)=soft_rx_bits;
       end
       
        % soft_rx_bits =[
        % 0
        % 0
        % 1
        % 1
        % 0
        % 1
        % 0
        % 1
        % 1 ]   
       soft_error_mapper = sum(abs(mimo_conv_encoded_bits - mimo_soft_rx_bits));
       soft_error_mappers(i_packet) = sum(soft_error_mapper);
       
       % re-organization bit structure change to properly format to send to
       
        % mimo_soft_rx_bits =
        % 0     1
        % 0     1
        % 1     1
        % 1     0
        % 0     0
        % 1     0
        % 0     0
        % 1     1
       mimo_viterbi_bits = zeros(2,N_frame*b);
       for rx_idx=1:row
           soft_rx_bits = mimo_soft_rx_bits(:,rx_idx);
           soft_rx_bits_reshape = reshape(soft_rx_bits,[2,length(soft_rx_bits)/2]);
           
            % viterbi_bits =
            % 0     1     0     1     1     1     0     1     0     1     0     0     1     1     0     0
           viterbi_bits = myviterbi(soft_rx_bits_reshape.');
           
           mimo_viterbi_bits(rx_idx,:) = viterbi_bits;
       end
       %keep error of each i packet, then we calculate BER
       soft_error = sum(abs(msg_symbol - mimo_viterbi_bits));
       if soft_error>0
           msg_symbol;
           mimo_viterbi_bits;
       end
       soft_errors(i_packet) = sum(soft_error);
    end %end for loop for i packet
    soft_mapper_BER(i_SNR) = sum(soft_error_mappers)/(N_packet*N_frame*b*2);
    soft_BER(i_SNR) = sum(soft_errors)/(N_packet*N_frame*b);
    
end %end for loop for i SNR


axes2 = axes('Parent',figure,'Yscale','log','YMinorTick','on','YminorGrid','on','FontSize',12,'FontName','Time New Roman');
xlim(axes2,[0 max(SNRdBs)]);
ylim(axes2,[1e-006 1]);
grid(axes2,'on');
hold(axes2,'on');
semilogy(SNRdBs, soft_mapper_BER, 'bo-')
hold on;
%Print SOFT decision BER vs SNR
semilogy(SNRdBs, soft_BER, 'ro-')
hold on;
grid on;
xlabel('SNR[dB]');
ylabel('BER');

clc
clear
close all
% begin
% N_packet = 1000; % No of iterations
N_packet = 3000;
b = 4; % modulation index 1:BPSK, 2:QPSK, 4: 16 QAM, 6: 64 QAM
N_frame = 1; % No of Modulation symbols per packet
M = 16;
SNRdBs = (10:5:40);
sq2 = sqrt(2);
mode = 'sdm';
channel = 'rayleigh';
plotconst = 'on';

for i_SNR = 1:length(SNRdBs)
    SNRdB = SNRdBs(i_SNR);
    %SNR=10^(SNRdB/10);
    SNR_minus=10^(-SNRdB/10);
    %sigma = sqrt(0.5/SNR);
    sigma = 10^(-SNRdB/20);
    for i_packet = 1:N_packet
       % Transmitter
       % build msg_symbol as bits block

       tx_error_symbols = [];
       hard_error_symbols = [];
       
       msg_symbol = round(rand([2 N_frame*b]));
       
       %% BEGIN CONVOLUTION ENCODING
       [row, col] = size(msg_symbol);
       mimo_conv_encoded_bits=zeros(col*2,row);
       
       for i=1:1:row
           conv_encoded_bits = convencode(msg_symbol(i,:));
           mimo_conv_encoded_bits(:,i)= conv_encoded_bits;
       end
       %% END CONVOLUTION ENCODING
       
       %% BEGIN 16QAM MAPPER AFTER CONVOLUTION ENCODING
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
       %% END 16QAM MAPPER AFTER CONVOLUTION ENCODING
        % mimo_tx_syms =
        %       X1                   X2
        % -3.0000 - 3.0000i  -3.0000 - 3.0000i
        % 1.0000 - 1.0000i   1.0000 + 3.0000i
       % Transmiter signal
       % X11 X12
       % X21 X22
       mimo_tx_syms = mimo_tx_syms.';
       
       X = mimo_tx_syms;
       
       %% BEGIN CHANNEL 
       H = zeros(N_frame*2,2); % RX antenal channel
       
       if strcmp(channel,'rayleigh')
           H = (randn(N_frame*2,2)+ 1i*randn(N_frame*2,2))/sq2;
       elseif strcmp(channel,'awgn')
           H = repmat([1,0],N_frame*2,2);
       else
           error('This channel is not supported');
       end
       %noise = signma*awgn
       N = (sigma*(randn(2,2)+1i*randn(2,2)))/sq2;
       %% END CHANNEL 
       
       %% RECEIVER
       Y =  H*X+N;
       
       Wzf = H^-1;
       mimo_rx_sym_zf = Wzf*Y;
       
       % Wmmse = (H'*H +1/SNR*I))^-1*H'
       Wmmse = (H'*H +(SNR_minus*eye(2)))^-1*H';
%        eye(2)/SNR
%        SNR_minus*eye(2)
       
       mimo_rx_sym_mmse = Wmmse*Y;
       
        % rx_sym =
        % -2.1500 - 3.6662i
        % -3.6876 + 0.4974i
        % 2.9912 - 5.3461i
        % -0.3573 + 3.2709i
        % -3.3309 + 2.1055i
        % 0.3879 - 1.0064i
        % 1.4661 - 1.5155i
        % -2.0477 + 4.0303i
       %hard_rx_bits = harddemapper(rx_sym);

        
        %% BEGIN ZERO FORCING CALCULATION EACH PACKET
        %%
       mimo_soft_rx_bits_zf = zeros(N_frame*b*2,2);
       
       for rx_idx=1:1:row
           rx_sym = mimo_rx_sym_zf(rx_idx,:);
           soft_rx_bits = softdemapper(rx_sym.');
           mimo_soft_rx_bits_zf(:,rx_idx)=soft_rx_bits;
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
       soft_error_mapper_zf = sum(abs(mimo_conv_encoded_bits - mimo_soft_rx_bits_zf));
       soft_error_mappers_zf(i_packet) = sum(soft_error_mapper_zf);
       
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
       mimo_viterbi_bits_zf = zeros(2,N_frame*b);
       for rx_idx=1:row
           soft_rx_bits = mimo_soft_rx_bits_zf(:,rx_idx);
           soft_rx_bits_reshape = reshape(soft_rx_bits,[2,length(soft_rx_bits)/2]);
           
            % viterbi_bits =
            % 0     1     0     1     1     1     0     1     0     1     0     0     1     1     0     0
           viterbi_bits = myviterbi(soft_rx_bits_reshape.');
           mimo_viterbi_bits_zf(rx_idx,:) = viterbi_bits;
       end
       %keep error of each i packet, then we calculate BER

       soft_error_zf = sum(abs(msg_symbol - mimo_viterbi_bits_zf));
       if sum(soft_error_zf)>0
           %X;
           %mimo_rx_sym;
           %msg_symbol;
           %mimo_viterbi_bits;
       end
       soft_errors_zf(i_packet) = sum(soft_error_zf);
       %% END ZERO FORCING CALCULATION EACH PACKET
       %%
       
       %% BEGIN MMSE CALCULATION EACH PACKET
       %%
       mimo_soft_rx_bits_mmse = zeros(N_frame*b*2,2);
       
       for rx_idx=1:1:row
           rx_sym = mimo_rx_sym_mmse(rx_idx,:);
           soft_rx_bits = softdemapper(rx_sym.');
           mimo_soft_rx_bits_mmse(:,rx_idx)=soft_rx_bits;
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
       soft_error_mapper_mmse = sum(abs(mimo_conv_encoded_bits - mimo_soft_rx_bits_mmse));
       soft_error_mappers_mmse(i_packet) = sum(soft_error_mapper_mmse);
       
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
       mimo_viterbi_bits_mmse = zeros(2,N_frame*b);
       for rx_idx=1:row
           soft_rx_bits = mimo_soft_rx_bits_mmse(:,rx_idx);
           soft_rx_bits_reshape = reshape(soft_rx_bits,[2,length(soft_rx_bits)/2]);
           
            % viterbi_bits =
            % 0     1     0     1     1     1     0     1     0     1     0     0     1     1     0     0
           viterbi_bits = myviterbi(soft_rx_bits_reshape.');
           mimo_viterbi_bits_mmse(rx_idx,:) = viterbi_bits;
       end
       %keep error of each i packet, then we calculate BER

       soft_error_mmse = sum(abs(msg_symbol - mimo_viterbi_bits_mmse));
       if sum(soft_error_mmse)>0
           %X;
           %mimo_rx_sym;
           %msg_symbol;
           %mimo_viterbi_bits;
       end
       soft_errors_mmse(i_packet) = sum(soft_error_mmse);
       %% END MMSE EACH CALCULATION PACKET
       %%
       
    end %end for loop for i packet
    
    % compute total error for zero forcing solution
    %soft_mapper_BER_zf(i_SNR) = sum(soft_error_mappers_zf)/(N_packet*N_frame*b*2);
    soft_BER_zf(i_SNR) = sum(soft_errors_zf)/(N_packet*N_frame*b);
    
    % compute total error for mmse solution
    %soft_mapper_BER_mmse(i_SNR) = sum(soft_error_mappers_mmse)/(N_packet*N_frame*b*2);
    soft_BER_mmse(i_SNR) = sum(soft_errors_mmse)/(N_packet*N_frame*b);
    
end %for loop for i SNR

axes2 = axes('Parent',figure,'Yscale','log','YMinorTick','on','YminorGrid','on','FontSize',12,'FontName','Time New Roman');
xlim(axes2,[0 max(SNRdBs)]);
ylim(axes2,[1e-006 1]);
grid(axes2,'on');
hold(axes2,'on');
% print ZERO FORCING BER
semilogy(SNRdBs, soft_BER_zf, 'bo-')
hold on;
% print MMSE BER
semilogy(SNRdBs, soft_BER_mmse, 'ro-')
hold on;
grid on;
xlabel('SNR[dB]');
ylabel('BER');
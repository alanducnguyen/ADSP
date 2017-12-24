clc
clear
close all
% begin
% N_packet = 1000; % No of iterations
N_packet = 10;
b = 4; % modulation index 1:BPSK, 2:QPSK, 4: 16 QAM, 6: 64 QAM
N_frame = 2; % No of Modulation symbols per packet
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

       mimo_tx_syms = zeros(N_frame,2);
       
       for idx=1:2
            % tx_sym =
            % -3.0000 + 1.0000i
            % 3.0000 + 1.0000i
            % 1.0000 - 1.0000i
            % -1.0000 - 3.0000i
            % -1.0000 + 3.0000i
            % -1.0000 - 3.0000i
            % 3.0000 + 1.0000i
            % -3.0000 + 1.0000i
            
           tx_sym = mapper(msg_symbol(idx,:),b,N_frame);
           mimo_tx_syms(:,idx) = tx_sym;
       end
       
        % mimo_tx_syms =
        %       X1                   X2
        % -3.0000 - 3.0000i  -3.0000 - 3.0000i
        % 1.0000 - 1.0000i   1.0000 + 3.0000i
       % Transmiter signal
       % X11 X12
       % X21 X22
       mimo_tx_syms;
       
       mimo_tx_syms = mimo_tx_syms.';
       
       X = mimo_tx_syms;
       
       X
       %channel
       H = zeros(N_frame,2); % RX antenal channel
      
       if strcmp(channel,'rayleigh')
           H = (randn(N_frame,2)+ 1i*randn(N_frame,2))/sq2;
       elseif strcmp(channel,'awgn')
           H = repmat([1,0],N_frame,2);
       else
           error('This channel is not supported');
       end
       
       %noise = signma*awgn
       N = sigma*(randn(N_frame,2)+randn(N_frame,2)*1i);
      
       %Receiver
       R =  H*X+N; 
     
       H_inv = H^-1;
       
       % X signal extraction after adding noise
       mimo_rx_sym = H_inv*R;
       
       mimo_soft_rx_bits = zeros(N_frame*b,2);
       
       for rx_idx=1:1:2
           rx_sym = mimo_rx_sym(rx_idx,:);
           soft_rx_bits = softdemapper(rx_sym.');
           mimo_soft_rx_bits(:,rx_idx)=soft_rx_bits;
       end
       
       soft_error_mapper = sum(abs(msg_symbol.' - mimo_soft_rx_bits));
       soft_error_mappers(i_packet) = sum(soft_error_mapper);
    end
    soft_BER(i_SNR) = sum(soft_error_mappers)/(N_packet*N_frame*b*2);
end

axes2 = axes('Parent',figure,'Yscale','log','YMinorTick','on','YminorGrid','on','FontSize',12,'FontName','Time New Roman');
xlim(axes2,[0 max(SNRdBs)]);
ylim(axes2,[1e-006 1]);
grid(axes2,'on');
hold(axes2,'on');
%Print SOFT decision BER vs SNR
semilogy(SNRdBs, soft_BER, 'ro-')
hold on;
grid on;
xlabel('SNR[dB]');
ylabel('BER');
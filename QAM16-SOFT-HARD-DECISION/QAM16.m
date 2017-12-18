clc
clear
close all
% begin
N_packet = 1000; % No of iterations
%N_packet = 2;
N_frame = 2; % No of Modulation symbols per packet
M = 16;
b = log2(M); % bits per symbol
SNRdBs = (0:5:40);

%Es_N0_dB  = Eb_N0_dB + 10*log10(k);
%SNRdBs  = SNRdBs + 10*log10(b);
channel = 'rayleigh';
plotconst = 'on';

%send N_packet with difference power,range +5.
for i_SNR = 1:length(SNRdBs)
    SNRdB = SNRdBs(i_SNR);
    sigma = sqrt(0.5/(10^(SNRdB/10)));
    for i_packet = 1:N_packet
       % Transmitter
       % build msg_symbol as bits block
       % 1     0     0     1
       % 0     0     1     0
       tx_error_symbols = [];
       hard_error_symbols = [];
       msg_symbol = randi([0 1],[N_frame*b,1]);
       tx_bits = msg_symbol;
       tx_sym = mapper(tx_bits.',b,N_frame);
     
       % Transmiter signal
       X = zeros(N_frame,1);
       X = tx_sym;
       %channel
       H = zeros(N_frame,1); % RX antenal channel
       if strcmp(channel,'rayleigh')
           H = (randn(N_frame,1)+ j*rand(N_frame,1))/sqrt(2);  % white guassian noise, 0dB variance
       elseif strcmp(channel,'awgn')
           H = repmat([1,0],N_frame,1);
       else
           error('This channel is not supported');
       end
       
       % H.*X is the element-by-element product of H and X.
       N = (randn(length(X),1)+randn(length(X),1)*1i)*sigma;
       R =  H.*X+N; %noise = signma*awgn
       %R =  H.*X; %without noise
       
       %Retrieve rx signal with noise
       rx_sym = R./H;
       
       hard_rx_bits = harddemapper(rx_sym);
       soft_rx_bits = softdemapper(rx_sym);

       %keep error of each i packet, then we calculate BER
       %hard_error = sum(abs(tx_bits - hard_rx_bits));
       
       hard_error = biterr(tx_bits,hard_rx_bits);
       if(hard_error>0)
           tx_error_symbols = [tx_error_symbols tx_sym];
           hard_error_symbols = [hard_error_symbols rx_sym];
       end
       soft_error = sum(abs(tx_bits - soft_rx_bits));
       
       hard_errors(i_packet) = sum(hard_error);
       soft_errors(i_packet) = sum(soft_error);
       
    end %end for loop for i packet
    
    tx_iphases = [];
    tx_qphases = [];
    
    rx_iphases = [];
    rx_qphases = [];
    has_error = length(tx_error_symbols);
    if has_error>=1
        figure;
        for idx=1:1:length(tx_error_symbols)
            tx_error_sym = tx_error_symbols(idx);
            rx_error_sym = hard_error_symbols(idx);
            if strcmp(plotconst ,'on')
                tx_iphase = real(tx_error_sym);
                tx_qphase = imag(tx_error_sym);
                rx_iphase = real(rx_error_sym);
                rx_qphase = imag(rx_error_sym);
                tx_iphases(idx) = tx_iphase;
                tx_qphases(idx) = tx_qphase;
                rx_iphases(idx) = rx_iphase;
                rx_qphases(idx) = rx_qphase;
            end
        end
        scatter(tx_iphases,tx_qphases,'*r');
        hold on
        scatter(rx_iphases,rx_qphases,'*b');
        axis([-5 5 -5 5]);
        grid on;
        hold off;
    end
    %hard_errors
    %sum(hard_errors)
    hard_BER(i_SNR) = sum(hard_errors)/(N_packet*N_frame*b);
    %hard_BER(i_SNR)
    soft_BER(i_SNR) = sum(soft_errors)/(N_packet*N_frame*b);
end %end for loop for i SNR

%theoryBer = 3/2*erfc(sqrt(0.1*(10.^(SNRdBs/10))));

figure;
axes2 = axes('Parent',figure,'Yscale','log','YMinorTick','on','YminorGrid','on','FontSize',12,'FontName','Time New Roman');
xlim(axes2,[0 max(SNRdBs)]);
ylim(axes2,[1e-006 1]);
grid(axes2,'on');
hold(axes2,'on');

%Print HARD decision BER vs SNR
%semilogy(SNRdBs, hard_BER,'Marker','+','LineWidth',1,'Color','r');
hold on;
%Print SOFT decision BER vs SNR
semilogy(SNRdBs, soft_BER, 'mx-','Linewidth',1);
%semilogy(SNRdBs, theoryBer, 'bo-','Linewidth',1);
grid on;
%axis([5 20 10^-5 5])
xlabel('SNR[dB]');
ylabel('BER');
%legend('hard demapper','soft demapper','theory');
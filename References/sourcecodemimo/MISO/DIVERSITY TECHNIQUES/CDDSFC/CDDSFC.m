
clc;
clear all;
close all;

N=1024;
j=sqrt(-1);
bpsk=randsrc(1,N);

%%----- Tx1,Tx2 Transmission strategy----- %% 
for i=1:2:N
    Tx1(i)=bpsk(i);
    Tx1(i+1)=-conj(bpsk(i+1));
    Tx2(i)=bpsk(i+1);
    Tx2(i+1)=conj(bpsk(i));
end

for Eb_No_dB=0:1:12
        sigamp1=sqrt(10.^(Eb_No_dB/10)).*Tx1;
        sigamp2=sqrt(10.^(Eb_No_dB/10)).*Tx2;
    for M=1:1:1
%%-----Channel-h1-----%%
        X1=randn(1,N/2);
        Y1=randn(1,N/2);
        unit_X1=X1/sqrt(var(X1));
        unit_Y1=Y1/sqrt(var(Y1));
        alpha1=sqrt((unit_X1.^2)+(unit_Y1.^2));
        alpha_norm1=alpha1/sqrt((var(alpha1)));
        phase1=rand(1,N/2);
        uni_phase1=2*pi*phase1;        
        e_ph1=complex(cos(uni_phase1),-sin(uni_phase1));        
        alpha_e_ph1=(alpha_norm1.*e_ph1)/sqrt(var(alpha_norm1.*e_ph1));
        h11=fft((alpha_e_ph1),N/2)/sqrt(N/2);%%---CHANNEL_1 FFT---%%
            for i=1:1:N/2
                h1((2*i)-1)=h11(i); %%Two consecutive subcarrier channel is constant%%
                h1(2*i)=h11(i);
            end

  %%-----Channel-h2------%%
        X2=randn(1,N/2);
        Y2=randn(1,N/2);
        unit_X2=X2/sqrt(var(X2));
        unit_Y2=Y2/sqrt(var(Y2));
        alpha2=sqrt((unit_X2.^2)+(unit_Y2.^2));
        alpha_norm2=alpha2/sqrt((mean(alpha2.^2)));
        phase2=rand(1,N/2);
        uni_phase2=2*pi*phase2;        
        e_ph2=complex(cos(uni_phase2),-sin(uni_phase2));
        alpha_e_ph2=(alpha_norm2.*e_ph2)/sqrt(var(alpha_norm2.*e_ph2));        
        h22=(exp(-j*2*pi*(1:N/2)*1/(N/2)).*fft((alpha_e_ph2),N/2))/sqrt(N/2);%%---CHANNEL_2 FFT---%%
        for i=1:1:N/2 
            h2((2*i)-1)=h22(i); %%Two consecutive subcarrier channel is constant%%       
            h2(2*i)=h22(i);            
        end

  %%-----AWGN------%%
        noise11=randn(1,N);
        noise12=randn(1,N);
        uni_var_noise11=noise11/sqrt(var(noise11));
        uni_var_noise12=noise12/sqrt(var(noise12));
        cgn=complex(uni_var_noise11,uni_var_noise12); %%--Varience of noise is 2 befor FFT--%%
        %var(cgn)
        cmplex_noise=fft(cgn,N)/sqrt(N);%%--Varience of noise is 2 after FFT--%%
        
        %var(cmplex_noise1)
        
 %%---------Receiver--------%%
 %%-----Received Signal-----%%
    for i=1:1:N/2
        r1((2*i)-1)=h1((2*i)-1)*sigamp1((2*i)-1)+h2((2*i)-1)*sigamp2((2*i)-1)+cmplex_noise((2*i)-1);    
        r1(2*i)=h1(2*i)*conj(sigamp1(2*i))+h2(2*i)*conj(sigamp2(2*i))+cmplex_noise(2*i);
    end
  %%---Signal Estimation from Received Signal---%%  
    for i=1:1:N/2
        sest((2*i)-1)=conj(h1((2*i)-1))*r1((2*i)-1)+h2(2*i)*conj(r1(2*i));
        sest(2*i)=conj(h2((2*i)-1))*r1((2*i)-1)-h1(2*i)*conj(r1(2*i));
    end
   %%---Estimated Signal Decoding---%%
        for i=1:1:N
             if ((dist(sest(i),-1))<(dist(sest(i),1)))
                decision(i)=-1;
            else
                decision(i)=1;
            end
        end
        %%---Bit Error Calculation---%
        error=0;        
        for i=1:1:N
                if (decision(i)~=bpsk(i))
                    error=error+1;
                end
        end
        Ber(M)=error/N;
    end
    Pb_ASFC_CDD(Eb_No_dB+1)=mean(Ber);
end        
Eb_No_dB=0:1:12;
semilogy(Eb_No_dB,Pb_ASFC_CDD,'r--');
grid on;
hold on;
axis([0 12 10^-4 0.5]);
legend('Alamouti Space Frequency Coding with Tx1_CD=0,Tx2_CD=1');
%legend('ASTBC','ASFBC');
xlabel('SNR(dB)');
ylabel('BER');
title('Alamouti Space Frequency Coding CDD (MT=2,MR=1)in Rayleigh fading channel');
%title('ASTBC,ASFBC (MT=2,MR=1)in Rayleigh fading channel')
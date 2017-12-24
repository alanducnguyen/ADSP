clc;
clear all;
close all;

N=1000;
bpsk=randsrc(1,N);
j=sqrt(-1);

%%----- Tx1,Tx2 Transmission strategy----- %% 

Tx=bpsk;

for Eb_No_dB=0:1:12
    sigamp=sqrt(10.^(Eb_No_dB/10)).*Tx;
  
    for M=1:1:10
%%-----Channel-h1-----%%
        X1=randn(1,N);
        Y1=randn(1,N);
        unit_X1=X1/sqrt(var(X1));
        unit_Y1=Y1/sqrt(var(Y1));
        alpha1=sqrt((unit_X1.^2)+(unit_Y1.^2));
        alpha_norm1=alpha1/sqrt((mean(alpha1.^2)));
        phase1=rand(1,N);
        uni_phase1=2*pi*phase1;
%         uni_phase_norm1=uni_phase1/sqrt(var(uni_phase1));
%         e_ph1=complex(cos(uni_phase_norm1),-sin(uni_phase_norm1));
        uni_phase_norm1=uni_phase1;
        e_ph1=complex(cos(uni_phase_norm1),-sin(uni_phase_norm1));
        h1=alpha_norm1.*e_ph1;
            
  %%-----Channel-h2------%%
        X2=randn(1,N);
        Y2=randn(1,N);
        unit_X2=X2/sqrt(var(X2));
        unit_Y2=Y2/sqrt(var(Y2));
        alpha2=sqrt((unit_X2.^2)+(unit_Y2.^2));
        alpha_norm2=alpha2/sqrt((mean(alpha2.^2)));
        phase2=rand(1,N);
        uni_phase2=2*pi*phase2;
%         uni_phase_norm2=uni_phase2/sqrt(var(uni_phase2));
%         e_ph2=complex(cos(uni_phase_norm2),-sin(uni_phase_norm2));
        uni_phase_norm2=uni_phase2;
        e_ph2=complex(cos(uni_phase_norm2),-sin(uni_phase_norm2));
        h2=alpha_norm2.*e_ph2;
 %%------AWGN------%%
        noise1=randn(1,N);
        noise2=randn(1,N);
        uni_var_noise1=noise1/sqrt(var(noise1));
        uni_var_noise2=noise2/sqrt(var(noise2));
        cmplex_noise=complex(uni_var_noise1,uni_var_noise2);
        CGN=cmplex_noise/sqrt(2);

 %%------TRANSMIT ANTENNA SELECTION USING MAXIMUM POWER-----%%                            
            p1=(h1.*conj(h1));
            p2=(h2.*conj(h2));    
            Rx1=h1.*sigamp;
            Rx2=h2.*sigamp;
                for i=1:1:N
                    if p1(i)==max(p1(i),p2(i))
                       h(i)=h1(i); 
                        Rx(i)=Rx1(i);  %%---Transmit Antenna1---%% 
                    else
                       h(i)=h2(i);   
                        Rx(i)=Rx2(i);  %%---Transmit Antenna2---%%
                    end
                end
 %%---------------------Receiver----------------------%%
                   
        Y=Rx+CGN;
        sest=conj(h).*Y;
        
            for i=1:1:N
                 if ((dist(sest(i),-1))<(dist(sest(i),1)))
                     decision(i)=-1;
                 else
                     decision(i)=1;
                 end
            end           
        
            error=0;        
            for i=1:1:N
                if (decision(i)~=bpsk(i))
                    error=error+1;
                end
            end
                    Ber(M)=error/N;
    end
    Pb_SDT(Eb_No_dB+1)=mean(Ber);
end

Eb_No_dB=0:1:12;
semilogy(Eb_No_dB,Pb_SDT,'r--');
grid on;
hold on;
axis([0 12 10^-4 0.5])
legend('SDT');
xlabel('SNR(dB)');
ylabel('BER');
title('Selection Diversity Transmission(MT=2,MR=1)in Rayleigh fading channel');

%end
        
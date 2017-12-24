% clc;
% clear all;
% close all;

function Pb_MRT=MRT(N)
%N=1000;
bpsk=randsrc(1,N);

%%----- Tx1,Tx2 Transmission strategy----- %% 

Tx1=bpsk;
Tx2=bpsk;

for Eb_No_dB=0:1:12
    sigamp1=sqrt(10.^(Eb_No_dB/10)).*Tx1;
    sigamp2=sqrt(10.^(Eb_No_dB/10)).*Tx2;
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
        uni_phase_norm1=uni_phase1/sqrt(var(uni_phase1));
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
        uni_phase_norm2=uni_phase2/sqrt(var(uni_phase2));
        e_ph2=complex(cos(uni_phase_norm2),-sin(uni_phase_norm2));
        h2=alpha_norm2.*e_ph2;
 %%------AWGN------%%
        noise11=randn(1,N);
        noise12=randn(1,N);
        uni_var_noise11=noise11/sqrt(var(noise11));
        uni_var_noise12=noise12/sqrt(var(noise12));
        cmplex_noise1=complex(uni_var_noise11,uni_var_noise12);
        CGN1=cmplex_noise1;
        noise21=randn(1,N);
        noise22=randn(1,N);
        uni_var_noise21=noise21/sqrt(var(noise21));
        uni_var_noise22=noise22/sqrt(var(noise22));
        cmplex_noise2=complex(uni_var_noise21,uni_var_noise22);
        CGN2=cmplex_noise2;
 %%------WEIGHT-----%%        
        
            w1=sqrt(2)*[conj(h1)./sqrt(abs(h1).^2+abs(h2).^2)];
            w2=sqrt(2)*[conj(h2)./sqrt(abs(h1).^2+abs(h2).^2)];
        
%         w1=w11./sqrt(var(w11));
%         w2=w22./sqrt(var(w22));
 %%---------Receiver--------%%
   
        Rx1=h1.*w1.*sigamp1;
        Rx2=h2.*w2.*sigamp2;
        
        Rx=Rx1+Rx2+CGN1;
        sest=Rx;
        
                       
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
    Pb_MRT(Eb_No_dB+1)=mean(Ber);
end
%MRT_F=MAT_MR(Pb_MRT);
% Eb_No_dB=0:1:12;
% semilogy(Eb_No_dB,Pb,'r--');
% grid on;
% hold on;
% axis([0 12 10^-4 0.5])
% legend('MRT');
% xlabel('SNR(dB)');
% ylabel('BER');
% title('Maximal Ratio Transmission(MT=2,MR=1)in Rayleigh fading channel');
end      








clc;
clear all;
close all;
%%--------Alamouti Scheme-Space Time Coding-MIMO-OFDM-------%%

N=1000;
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
    for M=1:1:1000
%%-----------Channel-h1----------%%

        X1=randn(1,N/2);
        Y1=randn(1,N/2);
        unit_X1=X1/sqrt(var(X1));
        unit_Y1=Y1/sqrt(var(Y1));
        alpha1=sqrt((unit_X1.^2)+(unit_Y1.^2));
        alpha_norm1=alpha1/sqrt((mean(alpha1.^2)));
        phase1=rand(1,N/2);
        uni_phase1=2*pi*phase1;           
        e_ph1=complex(cos(uni_phase1),-sin(uni_phase1));
        alpha_e_ph1=(alpha_norm1.*e_ph1)/sqrt(var(alpha_norm1.*e_ph1));        
        h11=fft((alpha_e_ph1),N/2)/sqrt(N/2);
        
            for i=1:1:N/2 
                h1((2*i)-1)=h11(i);
                h1(2*i)=h11(i);                 
            end

  %%------------Channel-h2-----------%%
  
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
        h22=fft((alpha_e_ph2),N/2)/sqrt(N/2);

        
        for i=1:1:N/2 
            h2((2*i)-1)=h22(i);            
            h2(2*i)=h22(i);
            
        end

   %%-----------Channel-h3------------%%

        X3=randn(1,N/2);
        Y3=randn(1,N/2);
        unit_X3=X3/sqrt(var(X3));
        unit_Y3=Y3/sqrt(var(Y3));
        alpha3=sqrt((unit_X3.^2)+(unit_Y3.^2));
        alpha_norm3=alpha3/sqrt((mean(alpha3.^2)));
        phase3=rand(1,N/2);
        uni_phase3=2*pi*phase3;        
        e_ph3=complex(cos(uni_phase3),-sin(uni_phase3));
        alpha_e_ph3=(alpha_norm3.*e_ph3)/sqrt(var(alpha_norm3.*e_ph3));        
        h33=fft((alpha_e_ph3),N/2)/sqrt(N/2);

        
            for i=1:1:N/2 
                h3((2*i)-1)=h33(i);
                h3(2*i)=h33(i);                 
            end
      
   %%-----------Channel-h4----------%%

        X4=randn(1,N/2);
        Y4=randn(1,N/2);
        unit_X4=X4/sqrt(var(X4));
        unit_Y4=Y4/sqrt(var(Y4));
        alpha4=sqrt((unit_X4.^2)+(unit_Y4.^2));
        alpha_norm4=alpha4/sqrt((mean(alpha4.^2)));
        phase4=rand(1,N/2);
        uni_phase4=2*pi*phase4;                
        e_ph4=complex(cos(uni_phase4),-sin(uni_phase4));
        alpha_e_ph4=(alpha_norm4.*e_ph4)/sqrt(var(alpha_norm1.*e_ph4));        
        h44=fft((alpha_e_ph4),N/2)/sqrt(N/2);

        
            for i=1:1:N/2 
                h4((2*i)-1)=h44(i);
                h4(2*i)=h44(i);                 
            end

  %%-------------AWGN--------------%%
        noise11=randn(1,N);
        noise12=randn(1,N);
        uni_var_noise11=noise11/sqrt(var(noise11));
        uni_var_noise12=noise12/sqrt(var(noise12));
        complex_noise1=complex(uni_var_noise11,uni_var_noise12);
        
        noise21=randn(1,N);
        noise22=randn(1,N);
        uni_var_noise21=noise21/sqrt(var(noise21));
        uni_var_noise22=noise22/sqrt(var(noise22));
        complex_noise2=complex(uni_var_noise21,uni_var_noise22);
        
 %%-------------Receiver-Rx1-----------%%
    for i=1:1:N
         r1(i)=h1(i)*sigamp1(i)+h2(i)*sigamp2(i)+complex_noise1(i);
    end
    for i=1:1:N/2
        sest1((2*i)-1)=(conj(h1((2*i)-1))*r1((2*i)-1))+(h2(2*i)*conj(r1(2*i)));
        sest1(2*i)=(conj(h2((2*i)-1))*r1((2*i)-1))-(h1(2*i)*conj(r1(2*i)));
    end
  %%-------------Receiver-Rx2-----------%%   
    for i=1:1:N
         r2(i)=h3(i)*sigamp1(i)+h4(i)*sigamp2(i)+complex_noise2(i);
    end
    for i=1:1:N/2
        sest2((2*i)-1)=(conj(h3((2*i)-1))*r2((2*i)-1))+(h4(2*i)*conj(r2(2*i)));
        sest2(2*i)=(conj(h4((2*i)-1))*r2((2*i)-1))-(h3(2*i)*conj(r2(2*i)));
    end
        sest=sest1+sest2;

    
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
    Pb(Eb_No_dB+1)=mean(Ber);
end
Eb_No_dB=0:1:12;
semilogy(Eb_No_dB,Pb,'r--');
grid on;
hold on;
axis([0 12 10^-6 0.5]);
legend('Alamouti scheme SFC (MT=2,MR=2)');
xlabel('SNR(dB)');
ylabel('BER');
title('Alamouti scheme SFC (MT=2,MR=2)in Rayleigh fading channel');
%         plot(abs(h1),'r');
%         hold on;
%         plot(abs(h2),'k');
%         hold on;
%         plot(abs(h3),'g');
%         hold on;
%         plot(abs(h4),'y');
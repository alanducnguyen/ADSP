channel = "rayleigh";
N_frame=4;
sq2 = sqrt(2);
SNR = 30;

SRNratio = SNR/10*log(10);
X = rand(1,2) + rand(1,2)*1i;

X = X';

sigma = sqrt(0.5/(10^(SRNratio/10)));

N = (randn(2,1)+randn(2,1)*1i)*sigma

H = (randn(2) + randn(2)*1i)/sqrt(2)

Y = H*X + N % with noise

Ynonoise = H*X ; %without noise

H_inv = H^-1;
X
X_nonoise = H_inv*Ynonoise - H_inv*N
X_withnoise = H_inv*Y - H_inv*N



%X_wn

channel = "rayleigh";
N_frame=4;
sq2 = sqrt(2);
SNRdB = 30;

X = rand(1,4) + rand(1,4)*1i;
X = X';
sigma = sqrt(0.5/(10^(SNRdB/10)));
N = (randn(4,1)+randn(4,1)*1i)*sigma

H = (randn(4) + randn(4)*1i)/sqrt(2)

Y = H*X + N
Yn0 = H*X ; %without noise

H_inv = inv(H);
X_wn = H_inv*Yn0 - H_inv*N;
X
X_wn

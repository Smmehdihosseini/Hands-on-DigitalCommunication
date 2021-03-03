% 2. Error Probability for BPSK Modulation with N Samples in Each Symbol
% Author: Seyed Mohammad Mehdi Hosseini (Smmehdihosseini@gmail.com)

clear;close all;clc;

num_bits=1000;
M=2; 
db_snr=-15:15;
N_Oversampling=[1 2 4 8 16];
ber_snr=zeros(1,length(db_snr));
for i_N=1:length(N_Oversampling)
for iteration_snr=1:length(db_snr)
err_totoal=0; iteration=0;
while iteration<500 
iteration=iteration+1;
bi=randi([0 1],1,num_bits);
sym=2*bi-1;
sym2=repmat(sym,N_Oversampling(i_N),1);
o_sym=reshape(sym2,[1,N_Oversampling(i_N)*length(sym)]);
E_sym=(o_sym* o_sym')/length(o_sym);
bpsk_sym=o_sym./(E_sym);

n0=10^(-(db_snr(iteration_snr)/10));
noise=sqrt(n0)*randn(1,N_Oversampling(i_N)*length(sym));
recx1=bpsk_sym+noise;
recx11=reshape(recx1,[N_Oversampling(i_N),length(sym)]);
o_recx1=sum(recx11,1);
bi1_hat=o_recx1>0;
num_err1=sum(xor(bi,bi1_hat));
err_totoal=err_totoal+num_err1;
end
ber_snr(iteration_snr)=err_totoal/(iteration*num_bits); 
end
semilogy(db_snr,ber_snr);hold on;grid on;
end
theory = qfunc(sqrt(10.^(db_snr./10)));
semilogy(db_snr,theory,'b-*');hold on;
legend('N=1','N=2','N=4','N=8','N=16','Theory Based');
xlabel('SNR(dB)'); ylabel('BER');
ylim([1e-4 1])
title('Bit Error Rate Versus SNR in AWGN Channel with Oversampling');
% 1. Error Probability for BPSK Modulation in AWGN  - Part4
% Author: Seyed Mohammad Mehdi Hosseini (Smmehdihosseini@gmail.com)

;close all;clc;


num_bits=1000;
M=2;
db_snr=0:12;
ber1_snr=zeros(1,length(db_snr));
for iteration_snr=1:length(db_snr)
bits=randi([0 1],1,num_bits);
symbol1=2*bits-1;
Er_sym1=(symbol1* symbol1')/length(symbol1);
bpsk_symbols1=symbol1./Er_sym1;

n0=10^(-(db_snr(iteration_snr)/10));
noise=sqrt(n0)*randn(1,length(symbol1));
recx1=bpsk_symbols1+noise;

if iteration_snr==3
scatterplot(bpsk_symbols1);
title('SNR=3 BPSK Symbols');
scatterplot(recx1);
title('SNR=3 Recieved');
elseif iteration_snr==8
scatterplot(bpsk_symbols1);
title('SNR=8 BPSK Symbols');
scatterplot(recx1);
title('SNR=8 Recieved');
elseif iteration_snr==13
scatterplot(bpsk_symbols1);
title('SNR=13 BPSK Symbols');
scatterplot(recx1);
title('SNR=13 Recieved');
end

bi1_hat=recx1>0;
ber1_snr(iteration_snr)=sum(xor(bits,bi1_hat))/num_bits;
end
% 1. Error Probability for BPSK Modulation in AWGN Channel - Part1
% Author: Seyed Mohammad Mehdi Hosseini (Smmehdihosseini@gmail.com)

clear;close all;clc;


num_bits=1000;
M=2;
db_snr=0:12;
ber1_snr=zeros(1,length(db_snr));
for iteration_snr=1:length(db_snr)
erroro_total1=0; iteration=0;
while iteration<100 || erroro_total1<100
iteration=iteration+1;
bits=randi([0 1],1,num_bits);
symbol1=2*bits-1;
Er_sym1=(symbol1* symbol1')/length(symbol1);
bpsk_symbols1=symbol1./Er_sym1;

n0=10^(-(db_snr(iteration_snr)/10));
noise=sqrt(n0)*randn(1,length(symbol1));
recx1=bpsk_symbols1+noise;




bi1_hat=real(recx1)>0;
num_error1=sum(xor(bits,bi1_hat));
erroro_total1=erroro_total1+num_error1;
end
ber1_snr(iteration_snr)=erroro_total1/(iteration*num_bits);
end
semilogy(db_snr,ber1_snr,'b-*');hold on;grid on;
theoryBer = qfunc(sqrt(10.^(db_snr./10)));
semilogy(db_snr,theoryBer,'g-');hold on;
title('Section 1 - Normal BPSK in AWGN');
legend('Simulation Value','Theory Value');
xlabel('SNR(dB)');
ylabel('BER');




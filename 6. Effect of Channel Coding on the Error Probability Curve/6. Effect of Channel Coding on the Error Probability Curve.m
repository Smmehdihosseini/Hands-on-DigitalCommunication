% 6. Effect of Channel Coding on the Error Probability Curve
% Author: Seyed Mohammad Mehdi Hosseini (Smmehdihosseini@gmail.com)



clear;close all;clc;

k=4;
n=7;
num_bits=1e6;
eb_n0_db=0:20;
es_n0_db_coding = eb_n0_db - 10*log10(n/k);
es_n0_db_noncoding = eb_n0_db;
n0_coding=10.^(-(es_n0_db_coding/10)); n0=10.^(-(es_n0_db_noncoding/10));
M=2;
P=[1 0 1;1 1 1;1 1 0;0 1 1];
G=[eye(k) P]; 
H=[P' eye(n-k)]; 
ber_cod=zeros(1,length(eb_n0_db)); ber=zeros(1,length(eb_n0_db));
for i_ebn0=1:length(eb_n0_db)
 b_in=randi([0 1],1,num_bits);
 bpsk = pskmod( b_in, M );
 noisedef = sqrt(n0(i_ebn0))*randn(1,length(bpsk));
 bpsk_noisy =( bpsk/sqrt((bpsk*bpsk')/length(bpsk)) )+ noisedef;
 bit_hat_pskmode = pskdemod(bpsk_noisy,M);
 blk_msg=zeros(num_bits/k,k);
 for message_count=1:num_bits/k
 blk_msg(message_count,:) = b_in(k*message_count-(k-1):k*message_count);
 end

 coded_message = mod(blk_msg*G,2);
 bpsk_codingword = pskmod( coded_message, M );
 bpsk_codingword_vector = reshape(bpsk_codingword',num_bits/k*n,1);
 noise_coding = sqrt(n0_coding(i_ebn0))*randn(length(bpsk_codingword_vector),1);

 bpsk_noisy_cod =(bpsk_codingword_vector/sqrt((bpsk_codingword_vector'*bpsk_codingword_vector)/length(bpsk_codingword_vector)) )+ noise_coding;
 bpsk_noisy_reshape_mat = (reshape(bpsk_noisy_cod,n,num_bits/k))';

 bit_hat_mat = pskdemod(bpsk_noisy_reshape_mat,M);
 syn= mod(bit_hat_mat*H',2); 
 codeword_hat= zeros(num_bits/k,n);
 for i_syn=1:num_bits/k
 if syn(i_syn,:)==zeros(1,n-k)
 codeword_hat(i_syn,:)= bit_hat_mat(i_syn,:);
 else
 a=0;
 while a==0
 for i=1:n
 comp = xor((syn(i_syn,:))',H(:,i));
 if comp==zeros(n-k,1)
 a=a+1;
 break
end
 end
 end
 err_position=i;
b=bit_hat_mat(i_syn,:);
 corr_codeword_hat=[b(1:err_position-1) not(b(1,err_position)) b(err_position+1:end)];
 codeword_hat(i_syn,:)=corr_codeword_hat;
 end
 end
 bit_hat_cod=(reshape(codeword_hat(:,1:k)',num_bits/k*k,1))';
 ber_cod(i_ebn0)=sum(xor(b_in,bit_hat_cod))/(num_bits);
 ber(i_ebn0)=sum(xor(b_in,bit_hat_pskmode))/(num_bits);
end

figure;
semilogy(eb_n0_db,ber_cod,'g-d','Linewidth',1.5);hold on;grid on;
semilogy(eb_n0_db,ber,'b-*','Linewidth',1.5);
legend('Hamming Coded','Not coded');
xlabel('Eb/N0(dB)');
ylabel('BER');
title('Error Vs SNR - BPSK Using Hamming Coding - C(7,4)')
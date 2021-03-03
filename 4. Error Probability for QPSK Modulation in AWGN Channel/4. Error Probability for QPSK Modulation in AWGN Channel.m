% 4. Error Probability for QPSK Modulation in AWGN Channel
% Author: Seyed Mohammad Mehdi Hosseini (Smmehdihosseini@gmail.com)


%Part 1

l_a=400000;
M=4;
a=randi([0,3],1,l_a);
SNR=-40:12 ;
len_SNR=length(SNR);

txSig=pskmod(a,M,pi/4);

p_err_vec=zeros(1,len_SNR);

for k=1:len_SNR
     rxSig=awgn(txSig,SNR(k));
     data_hat=pskdemod(rxSig,4,pi/4);
     n_err=sum(data_hat~=a);
     p_err=n_err/l_a;
     p_err_vec(1,k)=p_err;
end

plot(SNR,p_err_vec)
xlabel('SNR_{dB}')
ylabel('P_{error}')
title('Probability of QPSK error')


%Part 2

l_a=10000;
a=randi([0,3],1,l_a);
oversamp_vec=[1 2 4 8 16] ;
len_N=length(oversamp_vec);



SNR=-40:12 ; %dB
len_SNR=length(SNR);

p_err_mat=zeros(len_N,len_SNR);

for i=1:len_N
    a_star=repmat(a,oversamp_vec(i),1);	
    a_star=reshape(a_star,1,oversamp_vec(i)*l_a);	   
    txSig=pskmod(a_star,4,pi/4);	
    for k=1:len_SNR           	
        rxSig=awgn(txSig,SNR(k));	
        data_hat=reshape(rxSig,oversamp_vec(i),l_a);
        data_hatt=sum(data_hat,1);	     
        data_hat=pskdemod(data_hatt,4,pi/4);
        n_err=sum(data_hat~=a);
        p_err=n_err/l_a;
        p_err_mat(i,k)=p_err;
     end
end

plot(SNR,p_err_mat)
xlabel('SNR_{dB}')
ylabel('P_{error}')
title('Probability of error with oversampling')
legend('N=1','N=2','N=4','N=8','N=16')


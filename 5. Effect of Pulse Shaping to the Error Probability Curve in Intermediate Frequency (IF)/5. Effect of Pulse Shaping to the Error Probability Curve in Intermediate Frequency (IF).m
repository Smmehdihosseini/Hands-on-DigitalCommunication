% 5. Effect of Pulse Shaping to the Error Probability Curve in Intermediate Frequency (IF)
% Author: Seyed Mohammad Mehdi Hosseini (Smmehdihosseini@gmail.com)

clear;close all;clc;
num_bit=1000;
Time_bit=1e-3;
bit_rate=1/Time_bit;
M=2; 
k=log2(M); 
N=8;
Fs=N*(1/Time_bit);
db_snr=-10:2:10;
fc=3e3; 
fil_order=[127,255,1023]; 
BandWidth=[bit_rate/4,bit_rate,2*bit_rate]; 
wn=(2*BandWidth)/(Fs);
l_lp=fil_order+1; 
lp_delay=fil_order/2; 
lP_filter=[fir1(fil_order(3),wn(1));fir1(fil_order(3),wn(2));fir1(fil_order(3),wn(3))];
figure;freqz(lP_filter(1,:));
figure;freqz(lP_filter(2,:));
figure;freqz(lP_filter(3,:));
h1 = spectrum.welch('Hamming',100,5)
Hpsd1=psd(h1,lP_filter(1,:),'Fs',Fs)
h2 = spectrum.welch('Hamming',100,5)
Hpsd2=psd(h2,lP_filter(2,:),'Fs',Fs)
h3 = spectrum.welch('Hamming',100,5)
Hpsd3=psd(h3,lP_filter(3,:),'Fs',Fs)
figure;
p1=plot(Hpsd1);
set(p1,'LineWidth',1.5,'Color','y');hold on;
plt2=plot(Hpsd2);
set(plt2,'LineWidth',1.5,'Color','b');hold on;
plt3=plot(Hpsd3);
set(plt3,'LineWidth',1.5,'Color','r');
legend('BandWidth=0.25KHz','BandWidth=1KHz','BandWidth=2KHz');grid on;

span = 8;
sps = N;
beta=0.35;
h_rcosine = rcosdesign(beta,span,sps,'normal');
h_rcosine=h_rcosine/sqrt(sum(h_rcosine.^2));
h_rc_delay=span/2*sps; 
ber_snr_rc=zeros(length(BandWidth),length(db_snr));
for i_snr=1:length(db_snr)
err_total_rc=zeros(length(BandWidth),1); ite=0;
while ite<100 
bi=randi([0 1],1,num_bit);
sym1=2*bi-1;
bpsk_sym=upsample(sym1,N);
sym_rc1=conv(h_rcosine,bpsk_sym);
t=(0:length(sym_rc1)-1)/Fs;
sym_rc_pb=real(sym_rc1.*exp(1j*2*pi*fc*t));
bpsk_rc=sym_rc_pb/sqrt(sum(sym_rc_pb.^2)/length(sym_rc_pb));
n0=(10^(-(db_snr(i_snr)/10)));
noise=sqrt(n0)*randn(1,length(sym_rc1));
rx_rc =bpsk_rc +noise;
sym_rc_bb=real(rx_rc.*exp(-1j*2*pi*fc*t));
LPF_sym_rc_bb1=conv2(sym_rc_bb,lP_filter);
LPF_sym_rc_bb=LPF_sym_rc_bb1(:,floor(lp_delay(3))+1:end);

if i_snr==6 && ite==1
h9 = spectrum.welch('Hamming',100,5)
Hpsd9=psd(h9,sym_rc_pb,'Fs',Fs)
figure;plt9=plot(Hpsd9);
set(plt9,'LineWidth',1.5,'Color','b');
legend('Passband RC Shaped signal');
h4 = spectrum.welch('Hamming',140,60)
Hpsd4=psd(h4,rx_rc,'Fs',Fs)
h5 = spectrum.welch('Hamming',120,20)
Hpsd5=psd(h5,sym_rc_bb,'Fs',Fs)
h6 = spectrum.welch('Hamming',100,5)
Hpsd6=psd(h6,LPF_sym_rc_bb1(1,:),'Fs',Fs)
h7 = spectrum.welch('Hamming',100,5)
Hpsd7=psd(h7,LPF_sym_rc_bb1(2,:),'Fs',Fs)
h8 = spectrum.welch('Hamming',100,5)
Hpsd8=psd(h8,LPF_sym_rc_bb1(3,:),'Fs',Fs)
figure;plt4=plot(Hpsd4);
set(plt4,'LineWidth',1.5,'Color','r');
legend('Noisy Signal,Passband');
figure;
plt5=plot(Hpsd5);
set(plt5,'LineWidth',1.5,'Color','b');
legend('Noisy Signal,Baseband');
figure;
plt6=plot(Hpsd6);
set(plt6,'LineWidth',1.5,'Color','m');hold on;
plt7=plot(Hpsd7);
set(plt7,'LineWidth',1.5,'Color','g');hold on;
plt8=plot(Hpsd8);
set(plt8,'LineWidth',1.5,'Color','k');
legend('BandWidth=0.25KHz','BandWidth=1KHz','BandWidth=2KHz');
title('LP-Filtered Signal')
end
rx1_rc=conv2(h_rcosine,LPF_sym_rc_bb);
rx2_rc=rx1_rc( :,h_rc_delay*2+1:N:(num_bit*N)+(h_rc_delay*2+1)-1);
bi_hat_rc=rx2_rc>0;
num_err_rc=sum([xor(bi,bi_hat_rc(1,:));xor(bi,bi_hat_rc(2,:));xor(bi,bi_hat_rc(3,:))],2);
err_total_rc=err_total_rc+num_err_rc;
ite=ite+1;
end

ber_snr_rc(:,i_snr)=err_total_rc/(ite*num_bit);
end
figure;semilogy(db_snr,ber_snr_rc(1,:),'b-o');hold on;grid on;
semilogy(db_snr,ber_snr_rc(2,:),'r-s');hold on;
semilogy(db_snr,ber_snr_rc(3,:),'g-*');
legend('BW=0.25kHz','BW=1kHz','BW=2kHz');
xlabel('SNR');
ylabel('BER');
title('BER vs. SNR,Rb=1 kb/s, N_F=1024')

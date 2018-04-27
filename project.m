% ECE 414 Project
% By: Henrietta Odiete 
% Student number: 20482444
% Instructor: Professor Guang Gong 
 
clc 
clear
close all 
% Question 1a
% converts the bit sequence to a decimal number 
 
% A0 even and odd bits 
A0_evenbits=[0 0 0 1 1 1 0 1 0 0 0 1 0 0 1 0 1 1 1 0 0 0 1 0 1 1 1 0 1 1 0 1];
A0_oddbits= [0 1 0 0 1 0 0 0 0 1 0 0 0 1 1 1 1 0 1 1 0 1 1 1 1 0 1 1 1 0 0 0];
A0=zeros(1,32);
for b=1:32
   if(A0_evenbits(b)== 0 && A0_oddbits(b)==0)
       A0(b)=0;
   end
   
   if(A0_evenbits(b)== 0 && A0_oddbits(b)==1)
       A0(b)=1;
   end 

   if(A0_evenbits(b)== 1 && A0_oddbits(b)==0)
     A0(b)=2;
   end 
   
   if(A0_evenbits(b)== 1 && A0_oddbits(b)==1)
     A0(b)=3;
   end     
end 
% To generate other symbols A1-A10 circular left shift was applied to A0 
% used to compare output at reciver from the transmitter 
A1=circshift(A0,[1 -1]);
A2=circshift(A0,[1 -2]);
A3=circshift(A0,[1 -3]);
A4=circshift(A0,[1 -4]);
A5=circshift(A0,[1 -5]);
A6=circshift(A0,[1 -6]);
A7=circshift(A0,[1 -7]);
A8=circshift(A0,[1 -8]);
A9=circshift(A0,[1 -9]);
 
% modulates the sequence using 4QAM
Xk= zeros(1,32);
for f= 1:32 
  Xk(f)=  qammod(A0(f),4);
end 
 
%Ouestion 1b
% performs inverse fast fourier transform to get xn
xn= ifft(Xk);
 
 
% Question 1c
% adding cyclic prefix to the OFDM symbol xn 
% gets last eight symbols and adds it to the start
cyclic_xn=[xn(25:end) xn];
 
% Question 1d
% for other symbols cyclic prefix A1-A9
% circular shift of cyclic_xn to generate other OFDM symbols 
cyclic_xn1=circshift(cyclic_xn, [1 -1]);
cyclic_xn2=circshift(cyclic_xn, [1 -2]);
cyclic_xn3=circshift(cyclic_xn, [1 -3]);
cyclic_xn4=circshift(cyclic_xn, [1 -4]);
cyclic_xn5=circshift(cyclic_xn, [1 -5]);
cyclic_xn6=circshift(cyclic_xn, [1 -6]);
cyclic_xn7=circshift(cyclic_xn, [1 -7]);
cyclic_xn8=circshift(cyclic_xn, [1 -8]);
cyclic_xn9=circshift(cyclic_xn, [1 -9]);
 
% Question 1e
% Case 1: The time duration Ts is infinite
s=zeros(1,32);
s(1:32)=sum(Xk);
q=1:32;
% plots the response for infinte Ts 
% plot(q,s);
 
% Case 2:The time duration Ts is finite
% OFDM transmit time Ts= 4e-6
Xt = 0;
Ts = 4*10^(-6);  % Ts represents OFDM signal transmit time
t=0:100;
for r = 1:32:32
    Xt = Xt + Xk(r)*exp(2*1j*t*pi/Ts);
end
% plots the response for finite Ts
plot(t,Xt); 
 
 
 
% for X(f): spectrum of X(t) when Ts is not infinite 
% ?f=k/Ts
Xf=0;
for k = 1:32:32
    Xf = Xf + Xk(k)*sinc(Ts*(t.^-1-(k*k/Ts)));
end
%plots the frequency spectrum X(f)
figure();plot (t,Xf);
 
% impulse response of channel
% used to form zero forcing matrix
hn= [1 0.5 0.6 0.2 0.1 0 0 0 0 0 0 0 0 0 0];
% Question 3 
 
% Ouestion 3a- received signal yn
yn= conv(hn,cyclic_xn);
% Question 3b- spectrum of received signal Yk 
yk= fft(yn);
 
% Question 3c- Recover Xk 
% creates zero forcing matrix from the channel hn 
% for A0 symbol 
mn= tril(toeplitz(hn));
v= [0 0 0 0 0 0 0 1 0 0 0 0 0 0 0]';
bn= (mn\v)';
A0_recov=cconv(yn,bn,32);
A0_recov=fft(A0_recov);
% Question 3e- Demodulate Xk to get data 
A0_recov=qamdemod(A0_recov,4);
 
 
% Question 3d- Convert to serial stream  
% for A1 symbol 
yn1= conv(hn,cyclic_xn1);
yk1= fft(yn1);
A1_recov=cconv(yk1,bn,32);
A1_recov=qamdemod(A1_recov,4);
 
% for A2 symbol 
yn2= conv(hn,cyclic_xn2);
yk2= fft(yn2);
A2_recov=cconv(yk2,bn,32);
A2_recov=qamdemod(A2_recov,4);
 
% for A3 symbol
yn3= conv(hn,cyclic_xn3);
yk3= fft(yn3);
A3_recov=cconv(yk3,bn,32);
A3_recov=qamdemod(A3_recov,4);
 
% for A4 symbol
yn4= conv(hn,cyclic_xn4);
yk4= fft(yn4);
A4_recov=cconv(yk4,bn,32);
A4_recov=qamdemod(A4_recov,4);
 
 
% for A5 symbol
yn5= conv(hn,cyclic_xn5);
yk5= fft(yn5);
A5_recov=cconv(yk5,bn,32);
A5_recov=qamdemod(A5_recov,4);
 
% for A6 symbol
yn6= conv(hn,cyclic_xn6);
yk6= fft(yn6);
A6_recov=cconv(yk6,bn,32);
A6_recov=qamdemod(A6_recov,4);
 
% for A7 symbol
yn7= conv(hn,cyclic_xn7);
yk7= fft(yn7);
A7_recov=cconv(yk7,bn,32);
A7_recov=qamdemod(A7_recov,4);
 
% for A8 symbol
yn8= conv(hn,cyclic_xn8);
yk8= fft(yn8);
A8_recov=cconv(yk8,bn,32);
A8_recov=qamdemod(A8_recov,4);
 
% for A9 symbol
yn9= conv(hn,cyclic_xn9);
yk9= fft(yn9);
A9_recov=cconv(yk9,bn,32);
A9_recov=qamdemod(A9_recov,4);
 
% Question 3f Probablity of error 
number_symbols=320; % total number of decimal symbols sent 
symbol_error=zeros(1,10); %number of error for each sequence A0-A9 
 
for d=1:32
    if A0(d)~=A0_recov(d)
        symbol_error(1)=symbol_error(1)+1;
    end 
end 
 
for d=1:32
    if A1(d)~=A1_recov(d)
        symbol_error(2)=symbol_error(2)+1;
    end 
end 
for d=1:32
    if A2(d)~=A2_recov(d)
        symbol_error(3)=symbol_error(3)+1;
    end 
end 
 
for d=1:32
    if A3(d)~=A3_recov(d)
        symbol_error(4)=symbol_error(4)+1;
    end 
end 
 
for d=1:32
    if A4(d)~=A4_recov(d)
        symbol_error(5)=symbol_error(5)+1;
    end 
end 
 
for d=1:32
    if A5(d)~=A5_recov(d)
        symbol_error(6)=symbol_error(6)+1;
    end 
end 
 
for d=1:32
    if A6(d)~=A6_recov(d)
        symbol_error(7)=symbol_error(7)+1;
    end 
end 
 
for d=1:32
    if A7(d)~=A7_recov(d)
        symbol_error(8)=symbol_error(8)+1;
    end 
end 
 
for d=1:32
    if A8(d)~=A8_recov(d)
        symbol_error(9)=symbol_error(9)+1;
    end 
end
 
for d=1:32
    if A9(d)~=A9_recov(d)
        symbol_error(10)=symbol_error(10)+1;
    end 
end
 
P_e= 100*sum(symbol_error)/(number_symbols);
 


clc;
close all;
clear all;

Io=rgb2gray(imread("Lenna.png"));
%Io = im2double(Io);
info  = imfinfo("Lenna.png");
a = info.BitDepth;
[M N] = size(Io);

%ENCRYPTION ALGORITHM
%STEP 1
KR = zeros(1,M);
KC = zeros(1,N);
A  = 0:1:(2^a)-1;
n= (2^a)-1;
for i=1:M
    KR(i) = randperm(n,1);
end

for j=1:N
    KC(j) = randperm(n,1);
end
%STEP 2
ITERmax = 1;
ITER = 0;

%while(ITER~=ITERmax)
    
%STEP 3
ITER = ITER+1;
%STEP 4
alpha = zeros(1,N);
sum_alpha=0;

for i =1:M
    for j=1:N
        alpha(i) = alpha(i)+Io(i,j);
    end
end
M_alpha = zeros(1,M);
for i=1:M
    M_alpha(i) = mod(uint8(alpha(i)),2);
end
Iscr = zeros(M,N);
for i=1:M
    if M_alpha(i) == 0
        Io(i,:)=circshift(Io(i,:),KR(i));
        
    else
        Io(i,:)=circshift(Io(i,:),-1*KR(i));
    end
    %imshow(Io)
end

%STEP 5

beta = zeros(1,M);
sum_beta=0;

for j =1:N
    for i=1:M
        beta(i) = beta(i)+Io(i,j);
    end
end
M_beta = zeros(1,N);
for i=1:N
    M_beta(i) = mod(uint8(beta(i)),2);
end


for j=1:N
    if M_beta(j) == 0
        Io(:,j)=circshift(Io(:,j),KC(j));
    else
        Io(:,j)=circshift(Io(:,j),-1*KC(j));
    end
    %imshow(Io)
end


%STEP 6
%I1 = zeros(M,N);
%U_Io = uint8(Io);
for i=1:M/2
    for j=1:N
        I1(2*i-1,j) = bitxor(double(Io(2*i-1,j)),KC(j));
        KC = flip(KC);
        I1(2*i,j) = bitxor(double(Io(2*i,j)),KC(j));
    end
    %imshow((I1))
end

%IENC = zeros(M,N);
for i=1:M
    for j=1:N/2
        IENC(i,2*j-1) = bitxor((I1(i,2*j-1)),KR(j));
        KR = flip(KR);
        IENC(i,2*j) = bitxor((I1(i,2*j)),KR(j));
    end
    %imshow(IENC)
end


%DENCRYPTION ALGORITHM
%STEP 1

ITER  = 0;
%STEP 2
ITER = ITER + 1;
%STEP 3
for i=1:M
    for j=1:N/2
        I1(i,2*j-1) = bitxor((IENC(i,2*j-1)),KR(j));
        KR = flip(KR);
        I1(i,2*j) = bitxor((IENC(i,2*j)),KR(j));
    end
end
%STEP 4
for i=1:M/2
    for j=1:N
        Io(2*i-1,j) = bitxor(round(I1(2*i-1,j)),KC(j));
        KC = flip(KC);
        Io(2*i,j) = bitxor(round(I1(2*i,j)),KC(j));
    end
end

%STEP 5
for j =1:N
    for i=1:M
        beta(i) = beta(i)+Io(i,j);
    end
end

for i=1:N
    M_beta(i) = mod(uint8(beta(i)),2);
end


for j=1:N
    if M_beta(j) == 0
        Io(:,j)=circshift(Io(:,j),-1*KC(j));
    else
        Io(:,j)=circshift(Io(:,j),KC(j));
    end
    imshow(Io)
end

%STEP 6 
for i =1:M
    for j=1:N
        alpha(i) = alpha(i)+Io(i,j);
    end
end

for i=1:M
    M_alpha(i) = mod(uint8(alpha(i)),2);
end


for i=1:M
    if M_alpha(i) == 0
        Io(i,:)=circshift(Io(i,:),-1*KR(i));
    else
        Io(i,:)=circshift(Io(i,:),KR(i));
    end
    imshow(Io)
end

%end  
 





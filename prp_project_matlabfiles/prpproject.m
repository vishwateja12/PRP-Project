clear;
close all;
clc;
e = 1;
for SNRdb = 1:10:100
    phi = 0;
    fs = 10000;
    Ts=(1/fs);
    A = 1;
    f0 = 0.1980;
    t = 0:Ts:10;
    x = A*cos(2*pi*f0*t + phi);% equuation
    SNRdb10 = SNRdb/10;
    SNR = 10^(SNRdb10);
    variance = (A^2)/(2*SNR);%from given relation
    W = sqrt(variance).*randn(size(x)); %Gaussian white noise W
    x = x + W; %Add the noise
    figure(1);
    subplot(5,2,e);
    plot(t,x);
    xlabel("t------->");
    ylabel("INPUT SIGNAL");
    title("SIGNAL FOR DIFFERENT VALUES  OF SNRdb ");
    e = e+1;
end
%{
in above we found A by changing different values of SNRdb and plotted the
signal output from the figures we can say that A increases the signal is
least effected as SNRdb is directly proportional to A. Now we do know that
threshold SNRdb is in between the 30-40db's.

%}


% probability of error in bsc(noise channel)
%% initialising to store all the errors from a-f
PERROR_3 = zeros(1,7);

%% itterations for n =15,k=10,p=0.015 

% initialise
n = 15;
k = 10;
p = 0.015;
A = 2^n;
K = 2^k;
% no.of iterations for C
NN = 5;
%initialisation for storing p(n,k,p,C)
PERROR_4 = ones(1,NN);
% loop for C to fine min of p(n,k,p,C)
for jj = 1:1:NN
    
    % selection of random C from 2^n
    x = randperm(A,K);
    x = x - ones(1,K);
    c = de2bi(x,n);
    % initialise E = 0
    E = 0;
    % no. of iterations to find error
    N = 2000;
    
    for i = 1:1:N
        
        % selecting c_ from C
        t = randperm(K,1);
        c_ = c(t,:);
        y_ = c_;
        % initialising hamming distance
        d = zeros(1,K);
        % probability random
        f = rand(1,n);
        
        % finding y_ for probability of flipping p
        % [0,0,0,0],[0.12,0,0.97,0.8]
        for j = 1:1:n
           if f(j) <= p
            if y_(j) == 1
                y_(j)=0;
            elseif y_(j) == 0
                y_(j) =1;
            end
           end
        end
        
        % finding hamming distance with each code in C [00][11][-1,-1][1,1]
        for ii = 1:1:K
            cc = c(ii,:);
            len = abs(y_-cc);
            d(1,ii) = sum(len);
        end
        
        % finding minimum length
        [dmin,index] = min(d);
        % estimating c_cap
        c_cap = c(index,:);
        
        %finding hamming distance for c_ and c_cap
        diff = abs(c_cap-c_);
        b = sum(diff);
        
        % Indicator
        if b == 0
            I = 0;
        elseif b ~= 0
            I =1;
        end
        
        E = E +I;
        
    end
    
    % stores the values of p(n,k,p,C) for different C
    PERROR_4(1,jj) = E/N;
end
% p error in min condition of p(n,k,p,C)
CC = 1 - (p*log2(1/p) + (1-p)*log2(1/(1-p)));
rate = k/n;
disp('CHANNEL CAPACITY FOR n = 15,k=10,p=0.015');
disp(CC);
disp('RATE OF THE PROCESS IS');
disp(rate);
PERROR_3(1,1) = min(PERROR_4);

%% itterations for n =15,k=10,p=0.1 

% initialise
n = 15;
k = 10;
p = 0.1;
A = 2^n;
K = 2^k;
% no.of iterations for C
NN = 5;
%initialisation for storing p(n,k,p,C)
PERROR_4 = ones(1,NN);
% loop for C to fine min of p(n,k,p,C)
for jj = 1:1:NN
    
    % selection of random C from 2^n
    x = randperm(A,K);
    x = x - ones(1,K);
    c = de2bi(x,n);
    % initialise E = 0
    E = 0;
    % no. of iterations to find error
    N = 2000;
    
    for i = 1:1:N
        
        % selecting c_ from C
        t = randperm(K,1);
        c_ = c(t,:);
        y_ = c_;
        % initialising hamming distance
        d = zeros(1,K);
        % probability random
        f = rand(1,n);
        
        % finding y_ for probability of flipping p
        for j = 1:1:n
           if f(j) <= p
            if y_(j) == 1
                y_(j)=0;
            elseif y_(j) == 0
                y_(j) =1;
            end
           end
        end
        
        % finding hamming distance with each code in C
        for ii = 1:1:K
            cc = c(ii,:);
            len = abs(y_-cc);
            d(1,ii) = sum(len);
        end
        
        % finding minimum length
        [dmin,index] = min(d);
        % estimating c_cap
        c_cap = c(index,:);
        
        %finding hamming distance for c_ and c_cap
        diff = abs(c_cap-c_);
        b = sum(diff);
        
        % Indicator
        if b == 0
            I = 0;
        elseif b ~= 0
            I =1;
        end
        
        E = E +I;
        
    end
    
    % stores the values of p(n,k,p,C) for different C
    PERROR_4(1,jj) = E/N;
end
% p error in min condition of p(n,k,p,C)
CC = 1 - (p*log2(1/p) + (1-p)*log2(1/(1-p)));
rate = k/n;
disp('CHANNEL CAPACITY FOR n = 15,k=10,p=0.1');
disp(CC);
disp('RATE OF THE PROCESS IS');
disp(rate);
PERROR_3(1,2) = min(PERROR_4);
%% itterations for n =15,k=10,p=0.45 

% initialise
n = 15;
k = 10;
p = 0.45;
A = 2^n;
K = 2^k;
% no.of iterations for C
NN = 5;
%initialisation for storing p(n,k,p,C)
PERROR_4 = ones(1,NN);
% loop for C to fine min of p(n,k,p,C)
for jj = 1:1:NN
    
    % selection of random C from 2^n
    x = randperm(A,K);
    x = x - ones(1,K);
    c = de2bi(x,n);
    % initialise E = 0
    E = 0;
    % no. of iterations to find error
    N = 2000;
    
    for i = 1:1:N
        
        % selecting c_ from C
        t = randperm(K,1);
        c_ = c(t,:);
        y_ = c_;
        % initialising hamming distance
        d = zeros(1,K);
        % probability random
        f = rand(1,n);
        
        % finding y_ for probability of flipping p
        for j = 1:1:n
           if f(j) <= p
            if y_(j) == 1
                y_(j)=0;
            elseif y_(j) == 0
                y_(j) =1;
            end
           end
        end
        
        % finding hamming distance with each code in C
        for ii = 1:1:K
            cc = c(ii,:);
            len = abs(y_-cc);
            d(1,ii) = sum(len);
        end
        
        % finding minimum length
        [dmin,index] = min(d);
        % estimating c_cap
        c_cap = c(index,:);
        
        %finding hamming distance for c_ and c_cap
        diff = abs(c_cap-c_);
        b = sum(diff);
        
        % Indicator
        if b == 0
            I = 0;
        elseif b ~= 0
            I =1;
        end
        
        E = E +I;
        
    end
    
    % stores the values of p(n,k,p,C) for different C
    PERROR_4(1,jj) = E/N;
end
% p error in min condition of p(n,k,p,C)
CC = 1 - (p*log2(1/p) + (1-p)*log2(1/(1-p)));
rate = k/n;
disp('CHANNEL CAPACITY FOR n = 15,k=10,p=0.45');
disp(CC);
disp('RATE OF THE PROCESS IS');
disp(rate);
PERROR_3(1,3) = min(PERROR_4);
%% itterations for n =20,k=10,p=0.015 

% initialise
n = 20;
k = 10;
p = 0.015;
A = 2^n;
K = 2^k;
% no.of iterations for C
NN = 5;
%initialisation for storing p(n,k,p,C)
PERROR_4 = ones(1,NN);
% loop for C to fine min of p(n,k,p,C)
for jj = 1:1:NN
    
    % selection of random C from 2^n
    x = randperm(A,K);
    x = x - ones(1,K);
    c = de2bi(x,n);
    % initialise E = 0
    E = 0;
    % no. of iterations to find error
    N = 2000;
    
    for i = 1:1:N
        
        % selecting c_ from C
        t = randperm(K,1);
        c_ = c(t,:);
        y_ = c_;
        % initialising hamming distance
        d = zeros(1,K);
        % probability random
        f = rand(1,n);
        
        % finding y_ for probability of flipping p
        for j = 1:1:n
           if f(j) <= p
            if y_(j) == 1
                y_(j)=0;
            elseif y_(j) == 0
                y_(j) =1;
            end
           end
        end
        
        % finding hamming distance with each code in C
        for ii = 1:1:K
            cc = c(ii,:);
            len = abs(y_-cc);
            d(1,ii) = sum(len);
        end
        
        % finding minimum length
        [dmin,index] = min(d);
        % estimating c_cap
        c_cap = c(index,:);
        
        %finding hamming distance for c_ and c_cap
        diff = abs(c_cap-c_);
        b = sum(diff);
        
        % Indicator
        if b == 0
            I = 0;
        elseif b ~= 0
            I =1;
        end
        
        E = E +I;
        
    end
    
    % stores the values of p(n,k,p,C) for different C
    PERROR_4(1,jj) = E/N;
end
% p error in min condition of p(n,k,p,C)
CC = 1 - (p*log2(1/p) + (1-p)*log2(1/(1-p)));
rate = k/n;
disp('CHANNEL CAPACITY FOR n = 20,k=10,p=0.015');
disp(CC);
disp('RATE OF THE PROCESS IS');
disp(rate);
PERROR_3(1,4) = min(PERROR_4);
%% itterations for n =20,k=10,p=0.1 

% initialise
n = 20;
k = 10;
p = 0.1;
A = 2^n;
K = 2^k;
% no.of iterations for C
NN = 5;
%initialisation for storing p(n,k,p,C)
PERROR_4 = ones(1,NN);
% loop for C to fine min of p(n,k,p,C)
for jj = 1:1:NN
    
    % selection of random C from 2^n
    x = randperm(A,K);
    x = x - ones(1,K);
    c = de2bi(x,n);
    % initialise E = 0
    E = 0;
    % no. of iterations to find error
    N = 2000;
    
    for i = 1:1:N
        
        % selecting c_ from C
        t = randperm(K,1);
        c_ = c(t,:);
        y_ = c_;
        % initialising hamming distance
        d = zeros(1,K);
        % probability random
        f = rand(1,n);
        
        % finding y_ for probability of flipping p
        for j = 1:1:n
           if f(j) <= p
            if y_(j) == 1
                y_(j)=0;
            elseif y_(j) == 0
                y_(j) =1;
            end
           end
        end
        
        % finding hamming distance with each code in C
        for ii = 1:1:K
            cc = c(ii,:);
            len = abs(y_-cc);
            d(1,ii) = sum(len);
        end
        
        % finding minimum length
        [dmin,index] = min(d);
        % estimating c_cap
        c_cap = c(index,:);
        
        %finding hamming distance for c_ and c_cap
        diff = abs(c_cap-c_);
        b = sum(diff);
        
        % Indicator
        if b == 0
            I = 0;
        elseif b ~= 0
            I =1;
        end
        
        E = E +I;
        
    end
    
    % stores the values of p(n,k,p,C) for different C
    PERROR_4(1,jj) = E/N;
end
% p error in min condition of p(n,k,p,C)
CC = 1 - (p*log2(1/p) + (1-p)*log2(1/(1-p)));
rate = k/n;
disp('CHANNEL CAPACITY FOR n = 20,k=10,p=0.1');
disp(CC);
disp('RATE OF THE PROCESS IS');
disp(rate);
PERROR_3(1,5) = min(PERROR_4);
%% itterations for n =20,k=10,p=0.45 

% initialise
n = 20;
k = 10;
p = 0.45;
A = 2^n;
K = 2^k;
% no.of iterations for C
NN = 5;
%initialisation for storing p(n,k,p,C)
PERROR_4 = ones(1,NN);
% loop for C to fine min of p(n,k,p,C)
for jj = 1:1:NN
    
    % selection of random C from 2^n
    x = randperm(A,K);
    x = x - ones(1,K);
    c = de2bi(x,n);
    % initialise E = 0
    E = 0;
    % no. of iterations to find error
    N = 2000;
    
    for i = 1:1:N
        
        % selecting c_ from C
        t = randperm(K,1);
        c_ = c(t,:);
        y_ = c_;
        % initialising hamming distance
        d = zeros(1,K);
        % probability random
        f = rand(1,n);
        
        % finding y_ for probability of flipping p
        for j = 1:1:n
           if f(j) <= p
            if y_(j) == 1
                y_(j)=0;
            elseif y_(j) == 0
                y_(j) =1;
            end
           end
        end
        
        % finding hamming distance with each code in C
        for ii = 1:1:K
            cc = c(ii,:);
            len = abs(y_-cc);
            d(1,ii) = sum(len);
        end
        
        % finding minimum length
        [dmin,index] = min(d);
        % estimating c_cap
        c_cap = c(index,:);
        
        %finding hamming distance for c_ and c_cap
        diff = abs(c_cap-c_);
        b = sum(diff);
        
        % Indicator
        if b == 0
            I = 0;
        elseif b ~= 0
            I =1;
        end
        
        E = E +I;
        
    end
    
    % stores the values of p(n,k,p,C) for different C
    PERROR_4(1,jj) = E/N;
end
% p error in min condition of p(n,k,p,C)
CC = 1 - (p*log2(1/p) + (1-p)*log2(1/(1-p)));
rate = k/n;
disp('CHANNEL CAPACITY FOR n = 20,k=10,p=0.45');
disp(CC);
disp('RATE OF THE PROCESS IS');
disp(rate);
PERROR_3(1,6) = min(PERROR_4);
%% itterations for n =50,k=10,p=0.015 

% initialise
n = 50;
k = 10;
p = 0.015;
A = 2^n;
K = 2^k;
% no.of iterations for C
NN = 5;
%initialisation for storing p(n,k,p,C)
PERROR_4 = ones(1,NN);
% loop for C to fine min of p(n,k,p,C)
for jj = 1:1:NN
    
    % selection of random C from 2^n
    x = randperm(A,K);
    x = x - ones(1,K);
    c = de2bi(x,n);
    % initialise E = 0
    E = 0;
    % no. of iterations to find error
    N = 1000;
    
    for i = 1:1:N
        
        % selecting c_ from C
        t = randperm(K,1);
        c_ = c(t,:);
        y_ = c_;
        % initialising hamming distance
        d = zeros(1,K);
        % probability random
        f = rand(1,n);
        
        % finding y_ for probability of flipping p
        for j = 1:1:n
           if f(j) <= p
            if y_(j) == 1
                y_(j)=0;
            elseif y_(j) == 0
                y_(j) =1;
            end
           end
        end
        
        % finding hamming distance with each code in C
        for ii = 1:1:K
            cc = c(ii,:);
            len = abs(y_-cc);
            d(1,ii) = sum(len);
        end
        
        % finding minimum length
        [dmin,index] = min(d);
        % estimating c_cap
        c_cap = c(index,:);
        
        %finding hamming distance for c_ and c_cap
        diff = abs(c_cap-c_);
        b = sum(diff);
        
        % Indicator
        if b == 0
            I = 0;
        elseif b ~= 0
            I =1;
        end
        
        E = E +I;
        
    end
    
    % stores the values of p(n,k,p,C) for different C
    PERROR_4(1,jj) = E/N;
end
% p error in min condition of p(n,k,p,C)
CC = 1 - (p*log2(1/p) + (1-p)*log2(1/(1-p)));
rate = k/n;
disp('CHANNEL CAPACITY FOR n = 50,k=10,p=0.015');
disp(CC);
disp('RATE OF THE PROCESS IS');
disp(rate);
PERROR_3(1,7) = min(PERROR_4);
%% PLOTTING THE FIGURES OF PROBABILITY FROM a-f IN ORDER
figure(2);
r = 1:7;
stem(r,PERROR_3);
xlim([0,8]);
xlabel('inputs n,p,k');
ylabel('P_{ERROR}');
title('P_{ERROR} FOR DIFFERENT VALUES OF p,n,k');

%%

%{
In this section we try to create a hufmann code using probabilities of each
 symbol 
%}

%importing data
fileID = fopen('file1.txt','r');
s = fread(fileID);
fclose(fileID);


str = char(s);
punc = '!':'~';
string = [punc,' '];
long=sum(ismember(str,string));
r = 1;

% finding probabilities
for k=1:length(string)
  prob(k,1)=sum(ismember(str,string(k)))/long;
  if prob(k,1) ~= 0
      
      pstr(r) = string(k);
      pstr_prob(r) = prob(k,1); 
      r = r+1;
  end
end

% sorting in descending order
[prob_sorted,indices] = sort(pstr_prob,'descend');
pstr_sorted = pstr(indices);
pstr_sorted = double(pstr_sorted);
% initialising symbol probability vector
symprobvec =  [pstr_sorted;prob_sorted];
lensymprobvec = length(pstr_sorted);

% recursive function
codebook = huffmann(symprobvec',lensymprobvec);

% initialising code
symb(:,1) = char(codebook(:,1));

for ii = 1:1:length(symb)
    code{ii} = '';
    c = de2bi(codebook(ii,2),codebook(ii,3),'left-msb');
    str_c = num2str(c);
    str_c(isspace(str_c)) = '';
    code{ii} = strcat(str_c);
    
end

fileID = fopen('file1encode.txt','w');
for i = 1:1:length(str)
    for j = 1:1:length(symb)
        if str(i,1) == char(symb(j,1))
            c = de2bi(codebook(j,2),codebook(j,3),'left-msb');
            str_c = num2str(c);
            str_c(isspace(str_c)) = '';
            fprintf(fileID,'%s',str_c); 
            break;
        end
    end
end

fclose(fileID);

fileID = fopen('file1encode.txt','r');
bitstream = fscanf(fileID,'%s');
fclose(fileID);

msg = huffmandec(symb,code',bitstream);

fileID = fopen('file1decode.txt','w');
fprintf(fileID,'%c',msg);
fclose(fileID);

%{

HERE WE OBSERVE THE BINARY CHANNEL CODING ERROR WITH VARIOUS n,k,p VALUE
AND OBSERVING THE RESULTS OF PROBABILITY ERROR AND PROVIDING EXAMPLES FOR
SHANNON THEOREM

CHANNEL CAPACITY FOR n = 15,k=10,p=0.015
    0.8876

RATE OF THE PROCESS IS
    0.6667
RATE < CHANNEL CAPACITY IN THIS BINARY CHANNEL WE GET LOW PROBABILITY ERROR(p = 0.0535 )

CHANNEL CAPACITY FOR n = 15,k=10,p=0.1
    0.5310

RATE OF THE PROCESS IS
    0.6667
RATE > CHANNEL CAPACITY IN THIS BINARY CHANNEL WE GET HIGH PROBABILITY ERROR(p = 0.4565 )

CHANNEL CAPACITY FOR n = 15,k=10,p=0.45
    0.0072

RATE OF THE PROCESS IS
    0.6667
RATE > CHANNEL CAPACITY IN THIS BINARY CHANNEL WE GET HIGH PROBABILITY ERROR(p = 0.9955 )

CHANNEL CAPACITY FOR n = 20,k=10,p=0.015
    0.8876

RATE OF THE PROCESS IS
    0.5000
RATE < CHANNEL CAPACITY IN THIS BINARY CHANNEL WE GET LOW PROBABILITY(p = 0.0045 )

CHANNEL CAPACITY FOR n = 20,k=10,p=0.1
    0.5310

RATE OF THE PROCESS IS
    0.5000
RATE < CHANNEL CAPACITY IN THIS BINARY CHANNEL WE GET LOW PROBABILITY ERROR (p = 0.2410 )

CHANNEL CAPACITY FOR n = 20,k=10,p=0.45
    0.0072

RATE OF THE PROCESS IS
    0.5000
RATE > CHANNEL CAPACITY IN THIS BINARY CHANNEL WE GET HIGH PROBABILITY ERROR(p = 0.9970 )

CHANNEL CAPACITY FOR n = 50,k=10,p=0.015
    0.8876

RATE OF THE PROCESS IS
    0.2000
RATE << CHANNEL CAPACITY IN THIS BINARY CHANNEL WE GET LOW PROBABILITY ERROR(p = 0 )
%}
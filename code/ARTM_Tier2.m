clc,clear,close all;
%----------------------------------数据段-----------------------------------%
h = [4,5]/16;         %调制指数
L = 3;                %关联长度
M = 4;                %进制数
A = 1;                %幅度
data_number = 1000;   %产生符号数
sample_number = 16;     %每符号采样点数
R = 50;                 %符号率（波特率） 一秒钟内发送的符号数目
T = 1/R;                %符号周期
fs = R*sample_number;   %采样率 每秒从连续信号中提取并组成离散信号的采样个数
fc = 100;              %载波频率
shape = 'cos';        %成型波形（rec:矩形，cos:余弦）
SNR = 10*log10(10);             %信噪比
kk = 0;
data_ori = randi([0,M-1],1,data_number);     %产生M进制比特
data = 2*(data_ori+1)-1-M;                   %映射为双极性幅度

%---------------------------帧头参数（MSK信号）------------------------------%
hp = 0.5;
Ap = 1;
Mp = 2;
Lp = 1;
shape_p = 'rec';
L0 = 16;                     %DA序列长度
noise_len = 400;                   %帧前噪声长度
data_len = data_number+L0;   %一帧信号的长度
Nw = data_len*sample_number+noise_len;

%----------------M进制帧头MSK信号（M=2 L=1 h=0.5）的产生---------------------%
DA = (Mp-1)*ones(1,L0/4);
data_preamble = [-DA DA DA -DA];
[signal_preamble,~] = CPM_mod(hp,Lp,sample_number,length(data_preamble),T,Ap,shape_p,data_preamble);
S = signal_preamble(1:sample_number*L0);
figure %画基带波形
subplot(211)
plot(real(signal_preamble)),axis([0 L0*sample_number -1.2 1.2]);
xlabel('信号样本点'),ylabel('幅度/V')
title('传统调制法MSK信号同相基带波形')
subplot(212)
plot(imag(signal_preamble)),axis([0 L0*sample_number -1.2 1.2]);
xlabel('信号样本点'),ylabel('幅度/V')
title('传统调制法MSK信号正交基带波形')

%----------------------------PAM分解参数------------------------------------%
[Beta_k_i,D_k] = GetBetakiDk(L);
[g_kn,D_k] = PAM(h,L,M,sample_number,T);
D_k_new = zeros(2^(L-1)^2*3,2);
for i = 1:2^(L-1)^2*3
    for j = 1:2
        D_k_new(i,j) = length(g_kn{i,j})/sample_number;
    end
end
%[state_all,state_change,state_all_change1,state_all_change0] = state(h,L,M);
%-------------------------------产生基带信号---------------------------------%
[s_base,fi] = CPM_mod( h,L,sample_number,data_number,T,A,shape,data,fc); %基带信号，相位
t=T/sample_number:T/sample_number:T*(data_number+L-1);

figure
subplot(211)
plot(real(s_base)),axis([0 L0*sample_number -1.2 1.2]);
xlabel('信号样本点'),ylabel('幅度/V')
title('传统调制法ARTM Tier2信号同相基带波形')
subplot(212)
plot(imag(s_base)),axis([0 L0*sample_number -1.2 1.2]);
xlabel('信号样本点'),ylabel('幅度/V')
title('传统调制法ARTM Tier2信号正交基带波形')

s = [signal_preamble s_base];
t = 0:1/fs:(length(s)-1)/fs;

%s=A*cos(2*pi*fc*t+fi);                                                  %加载波
s = s.*exp(1j*2*pi*fc*t);

s1 = real(s);
s2 = awgn(s1,SNR);

[fi_i,fi_q,b] = BaseRecovery(fc,fs,t,s2);   %滤波
s_r = fi_i+1i*fi_q;%接收信号
delay = 31;        %延迟
for i = 1:length(s_r)-delay
    s_r(1,i) = s_r(1,delay+i-1);
end
mm = randn(1,noise_len); %帧头前噪声段
sigpower = sum(s1.*conj(s1))/length(s1);        %信号平均功率
EbNo = 10^(SNR/10);                                         %将SNR（信噪比）的dB值转换成数值
sigma2 = sigpower/EbNo;
noise1 = sqrt(sigma2/2).*mm;                                %零均值，方差为sigma2的实高斯白噪声
signal_crb = [noise1 s_r];
% figure,plot(real(signal_crb))
% xlabel('信号样本点'),ylabel('幅值/V'),title('一帧信号')

%---------------------------------帧检测-----------------------------------%
L0 = sample_number*L0;
r = zeros(1,length(signal_crb)-L0);
for k = 1:length(signal_crb)-L0
    C = Nw-k;
    temp = 0;
    for d = 1:L0-1
        temp = temp+abs(sum(conj(signal_crb(k:L0+k-d-1)).*signal_crb(k+d:L0+k-1).*S(1:L0-d).*conj(S(d:L0-1))));
    end
    r(k) = C*(sum(abs(signal_crb(k:Nw-1)).^2)+2*temp);
end
% figure
% plot(r),xlabel('接收信号的采样点'),ylabel('δ'),title('帧检测结果')
position = find(r==max(r));
s_r = signal_crb(position+L0+1:end); %帧同步

r_k = zeros(2^(L-1)^2*3,1);
h_dem = zeros(1,data_number);      %每个符号对应的调制指数矩阵
for i = 1:data_number
    h_dem(1,i) = h(mod(i-1,2)+1);
end

a_k = zeros(1,48);
delay = 0;%延迟过就不延迟了
%--------------------------viterbi译码-------------------------------------%
n_dem = 1;
for ii=1:2^(L-1)^2*3
    kk=delay+((n_dem-1)*sample_number+1:(n_dem+D_k_new(ii)-1)*sample_number);
    r_k(ii,1)=sum(s_r(kk).*g_kn{ii,1});
end

sum_index1 = 1*h(mod(n_dem-1,2)+1);
sum_index0 = 1*h(mod(n_dem-1,2)+1);
b1_0_n = exp(1i*pi*2*(sum_index1));
b1_0_n_1 = 1;
b1_0_n_2 = 1;
b1_0_n_3 = 1;
b1_1_n = exp(1i*pi*2*sum_index1);
b1_1_n_1 = 1;
b1_2_n = exp(1i*pi*2*sum_index1);
b1_3_n = exp(1i*pi*2*sum_index1);
b0_0_n = exp(1i*pi*(sum_index0));
b0_0_n_1 = 1;
b0_0_n_2 = 1;
b0_0_n_3 = 1;
b0_1_n = exp(1i*pi*sum_index0);
b0_1_n_1 = 1;
b0_2_n = exp(1i*pi*sum_index0);
b0_3_n = exp(1i*pi*sum_index0);

[a_k] = calculate( b1_0_n,b1_0_n_1,b1_0_n_2,b1_0_n_3,b1_1_n,b1_1_n_1,b1_2_n,b1_3_n,b0_0_n,b0_0_n_1,b0_0_n_2,b0_0_n_3,b0_1_n,b0_1_n_1,b0_2_n,b0_3_n);
value_new1=(real(conj(a_k)*r_k));  %传3路径值

sum_index1 = 1*h(mod(n_dem-1,2)+1);
sum_index0 = -1*h(mod(n_dem-1,2)+1);
b1_0_n = exp(1i*pi*2*(sum_index1));
b1_0_n_1 = 1;
b1_0_n_2 = 1;
b1_0_n_3 = 1;
b1_1_n = exp(1i*pi*2*sum_index1);
b1_1_n_1 = 1;
b1_2_n = exp(1i*pi*2*sum_index1);
b1_3_n = exp(1i*pi*2*sum_index1);
b0_0_n = exp(1i*pi*(sum_index0));
b0_0_n_1 = 1;
b0_0_n_2 = 1;
b0_0_n_3 = 1;
b0_1_n = exp(1i*pi*sum_index0);
b0_1_n_1 = 1;
b0_2_n = exp(1i*pi*sum_index0);
b0_3_n = exp(1i*pi*sum_index0);

[a_k] = calculate( b1_0_n,b1_0_n_1,b1_0_n_2,b1_0_n_3,b1_1_n,b1_1_n_1,b1_2_n,b1_3_n,b0_0_n,b0_0_n_1,b0_0_n_2,b0_0_n_3,b0_1_n,b0_1_n_1,b0_2_n,b0_3_n);
value_new2=(real(conj(a_k)*r_k));  %传1路径值

sum_index1 = -1*h(mod(n_dem-1,2)+1);
sum_index0 = 1*h(mod(n_dem-1,2)+1);
b1_0_n = exp(1i*pi*2*(sum_index1));
b1_0_n_1 = 1;
b1_0_n_2 = 1;
b1_0_n_3 = 1;
b1_1_n = exp(1i*pi*2*sum_index1);
b1_1_n_1 = 1;
b1_2_n = exp(1i*pi*2*sum_index1);
b1_3_n = exp(1i*pi*2*sum_index1);
b0_0_n = exp(1i*pi*(sum_index0));
b0_0_n_1 = 1;
b0_0_n_2 = 1;
b0_0_n_3 = 1;
b0_1_n = exp(1i*pi*sum_index0);
b0_1_n_1 = 1;
b0_2_n = exp(1i*pi*sum_index0);
b0_3_n = exp(1i*pi*sum_index0);

[a_k] = calculate( b1_0_n,b1_0_n_1,b1_0_n_2,b1_0_n_3,b1_1_n,b1_1_n_1,b1_2_n,b1_3_n,b0_0_n,b0_0_n_1,b0_0_n_2,b0_0_n_3,b0_1_n,b0_1_n_1,b0_2_n,b0_3_n);
value_new3=(real(conj(a_k)*r_k));  %传-1路径值

sum_index1 = -1*h(mod(n_dem-1,2)+1);
sum_index0 = -1*h(mod(n_dem-1,2)+1);
b1_0_n = exp(1i*pi*2*(sum_index1));
b1_0_n_1 = 1;
b1_0_n_2 = 1;
b1_0_n_3 = 1;
b1_1_n = exp(1i*pi*2*sum_index1);
b1_1_n_1 = 1;
b1_2_n = exp(1i*pi*2*sum_index1);
b1_3_n = exp(1i*pi*2*sum_index1);
b0_0_n = exp(1i*pi*(sum_index0));
b0_0_n_1 = 1;
b0_0_n_2 = 1;
b0_0_n_3 = 1;
b0_1_n = exp(1i*pi*sum_index0);
b0_1_n_1 = 1;
b0_2_n = exp(1i*pi*sum_index0);
b0_3_n = exp(1i*pi*sum_index0);


[a_k] = calculate( b1_0_n,b1_0_n_1,b1_0_n_2,b1_0_n_3,b1_1_n,b1_1_n_1,b1_2_n,b1_3_n,b0_0_n,b0_0_n_1,b0_0_n_2,b0_0_n_3,b0_1_n,b0_1_n_1,b0_2_n,b0_3_n);

value_new4=(real(conj(a_k)*r_k));  %传-3路径值

value1 = [value_new1;value_new2;value_new3;value_new4];%1时刻的路径值

k = 0;
for data_dem = 3:-2:-3
    k = k+1;
    if data_dem == 3     %多进制分解
        data_dem_new = [1,1];
    elseif data_dem == 1
        data_dem_new = [1,-1];
    elseif data_dem == -1
        data_dem_new = [-1,1];
    else
        data_dem_new = [-1,-1];
    end
    data_dem_true1 = data_dem_new(1,1);%gama_i,1
    data_dem_true0 = data_dem_new(1,2);%gama_i,0
    
    n_dem = 2;
    
    for ii=1:2^(L-1)^2*3
        kk=delay+((n_dem-1)*sample_number+1:(n_dem+D_k_new(ii)-1)*sample_number);
        r_k(ii,1)=sum(s_r(kk).*g_kn{ii,1});
    end
    
    sum_index1 = data_dem_true1*h(mod(n_dem,2)+1)+1*h(mod(n_dem-1,2)+1);
    sum_index0 = data_dem_true0*h(mod(n_dem,2)+1)+1*h(mod(n_dem-1,2)+1);
    b1_0_n = exp(1i*pi*2*(sum_index1));
    b1_0_n_1 = exp(1i*pi*2*(sum_index1-1*h(mod(n_dem-1,2)+1)));
    b1_0_n_2 = exp(1i*pi*2*(sum_index1-1*h(mod(n_dem-1,2)+1)-data_dem_true1*h(mod(n_dem,2)+1)));
    b1_0_n_3 = 1;
    b1_1_n = exp(1i*pi*2*(sum_index1-data_dem_true1*h(mod(n_dem,2)+1)));
    b1_1_n_1 = exp(1i*pi*2*(sum_index1-1*h(mod(n_dem-1,2)+1)));
    b1_2_n = exp(1i*pi*2*sum_index1);
    b1_3_n = exp(1i*pi*2*(sum_index1-data_dem_true1*h(mod(n_dem,2)+1)));
    b0_0_n = exp(1i*pi*(sum_index0));
    b0_0_n_1 = exp(1i*pi*(sum_index0)-1*h(mod(n_dem-1,2)+1));
    b0_0_n_2 = exp(1i*pi*(sum_index1-1*h(mod(n_dem-1,2)+1)-data_dem_true0*h(mod(n_dem,2)+1)));
    b0_0_n_3 = 1;
    b0_1_n = exp(1i*pi*sum_index0-data_dem_true0*h(mod(n_dem,2)+1));
    b0_1_n_1 = exp(1i*pi*(sum_index1-1*h(mod(n_dem-1,2)+1)));
    b0_2_n = exp(1i*pi*sum_index0);
    b0_3_n = exp(1i*pi*sum_index0-data_dem_true0*h(mod(n_dem,2)+1));
    
    [a_k] = calculate( b1_0_n,b1_0_n_1,b1_0_n_2,b1_0_n_3,b1_1_n,b1_1_n_1,b1_2_n,b1_3_n,b0_0_n,b0_0_n_1,b0_0_n_2,b0_0_n_3,b0_1_n,b0_1_n_1,b0_2_n,b0_3_n);
    value_new5=(real(conj(a_k)*r_k));  %传3路径值
    
    sum_index1 = data_dem_true1*h(mod(n_dem,2)+1)+1*h(mod(n_dem-1,2)+1);
    sum_index0 = data_dem_true0*h(mod(n_dem,2)+1)-1*h(mod(n_dem-1,2)+1);
    b1_0_n = exp(1i*pi*2*(sum_index1));
    b1_0_n_1 = exp(1i*pi*2*(sum_index1-1*h(mod(n_dem-1,2)+1)));
    b1_0_n_2 = exp(1i*pi*2*(sum_index1-1*h(mod(n_dem-1,2)+1)-data_dem_true1*h(mod(n_dem,2)+1)));
    b1_0_n_3 = 1;
    b1_1_n = exp(1i*pi*2*(sum_index1-data_dem_true1*h(mod(n_dem,2)+1)));
    b1_1_n_1 = exp(1i*pi*2*(sum_index1-1*h(mod(n_dem-1,2)+1)));
    b1_2_n = exp(1i*pi*2*sum_index1);
    b1_3_n = exp(1i*pi*2*(sum_index1-data_dem_true1*h(mod(n_dem,2)+1)));
    b0_0_n = exp(1i*pi*(sum_index0));
    b0_0_n_1 = exp(1i*pi*(sum_index0)+1*h(mod(n_dem-1,2)+1));
    b0_0_n_2 = exp(1i*pi*(sum_index1+1*h(mod(n_dem-1,2)+1)-data_dem_true0*h(mod(n_dem,2)+1)));
    b0_0_n_3 = 1;
    b0_1_n = exp(1i*pi*sum_index0-data_dem_true0*h(mod(n_dem,2)+1));
    b0_1_n_1 = exp(1i*pi*(sum_index1+1*h(mod(n_dem-1,2)+1)));
    b0_2_n = exp(1i*pi*sum_index0);
    b0_3_n = exp(1i*pi*sum_index0-data_dem_true0*h(mod(n_dem,2)+1));
    
    [a_k] = calculate( b1_0_n,b1_0_n_1,b1_0_n_2,b1_0_n_3,b1_1_n,b1_1_n_1,b1_2_n,b1_3_n,b0_0_n,b0_0_n_1,b0_0_n_2,b0_0_n_3,b0_1_n,b0_1_n_1,b0_2_n,b0_3_n);
    value_new6=(real(conj(a_k)*r_k));  %传1路径值
    
    sum_index1 = data_dem_true1*h(mod(n_dem,2)+1)-1*h(mod(n_dem-1,2)+1);
    sum_index0 = data_dem_true0*h(mod(n_dem,2)+1)+1*h(mod(n_dem-1,2)+1);
    b1_0_n = exp(1i*pi*2*(sum_index1));
    b1_0_n_1 = exp(1i*pi*2*(sum_index1+1*h(mod(n_dem-1,2)+1)));
    b1_0_n_2 = exp(1i*pi*2*(sum_index1+1*h(mod(n_dem-1,2)+1)-data_dem_true1*h(mod(n_dem,2)+1)));
    b1_0_n_3 = 1;
    b1_1_n = exp(1i*pi*2*(sum_index1-data_dem_true1*h(mod(n_dem,2)+1)));
    b1_1_n_1 = exp(1i*pi*2*(sum_index1+1*h(mod(n_dem-1,2)+1)));
    b1_2_n = exp(1i*pi*2*sum_index1);
    b1_3_n = exp(1i*pi*2*(sum_index1-data_dem_true1*h(mod(n_dem,2)+1)));
    b0_0_n = exp(1i*pi*(sum_index0));
    b0_0_n_1 = exp(1i*pi*(sum_index0)-1*h(mod(n_dem-1,2)+1));
    b0_0_n_2 = exp(1i*pi*(sum_index1-1*h(mod(n_dem-1,2)+1)-data_dem_true0*h(mod(n_dem,2)+1)));
    b0_0_n_3 = 1;
    b0_1_n = exp(1i*pi*sum_index0-data_dem_true0*h(mod(n_dem,2)+1));
    b0_1_n_1 = exp(1i*pi*(sum_index1-1*h(mod(n_dem-1,2)+1)));
    b0_2_n = exp(1i*pi*sum_index0);
    b0_3_n = exp(1i*pi*sum_index0-data_dem_true0*h(mod(n_dem,2)+1));
    
    [a_k] = calculate( b1_0_n,b1_0_n_1,b1_0_n_2,b1_0_n_3,b1_1_n,b1_1_n_1,b1_2_n,b1_3_n,b0_0_n,b0_0_n_1,b0_0_n_2,b0_0_n_3,b0_1_n,b0_1_n_1,b0_2_n,b0_3_n);
    value_new7=(real(conj(a_k)*r_k));  %传-1路径值
    
    sum_index1 = data_dem_true1*h(mod(n_dem,2)+1)-1*h(mod(n_dem-1,2)+1);
    sum_index0 = data_dem_true0*h(mod(n_dem,2)+1)-1*h(mod(n_dem-1,2)+1);
    b1_0_n = exp(1i*pi*2*(sum_index1));
    b1_0_n_1 = exp(1i*pi*2*(sum_index1+1*h(mod(n_dem-1,2)+1)));
    b1_0_n_2 = exp(1i*pi*2*(sum_index1+1*h(mod(n_dem-1,2)+1)-data_dem_true1*h(mod(n_dem,2)+1)));
    b1_0_n_3 = 1;
    b1_1_n = exp(1i*pi*2*(sum_index1-data_dem_true1*h(mod(n_dem,2)+1)));
    b1_1_n_1 = exp(1i*pi*2*(sum_index1+1*h(mod(n_dem-1,2)+1)));
    b1_2_n = exp(1i*pi*2*sum_index1);
    b1_3_n = exp(1i*pi*2*(sum_index1-data_dem_true1*h(mod(n_dem,2)+1)));
    b0_0_n = exp(1i*pi*(sum_index0));
    b0_0_n_1 = exp(1i*pi*(sum_index0)+1*h(mod(n_dem-1,2)+1));
    b0_0_n_2 = exp(1i*pi*(sum_index1+1*h(mod(n_dem-1,2)+1)-data_dem_true0*h(mod(n_dem,2)+1)));
    b0_0_n_3 = 1;
    b0_1_n = exp(1i*pi*sum_index0-data_dem_true0*h(mod(n_dem,2)+1));
    b0_1_n_1 = exp(1i*pi*(sum_index1+1*h(mod(n_dem-1,2)+1)));
    b0_2_n = exp(1i*pi*sum_index0);
    b0_3_n = exp(1i*pi*sum_index0-data_dem_true0*h(mod(n_dem,2)+1));
    
    [a_k] = calculate( b1_0_n,b1_0_n_1,b1_0_n_2,b1_0_n_3,b1_1_n,b1_1_n_1,b1_2_n,b1_3_n,b0_0_n,b0_0_n_1,b0_0_n_2,b0_0_n_3,b0_1_n,b0_1_n_1,b0_2_n,b0_3_n);
    
    value_new8=(real(conj(a_k)*r_k));  %传-3路径值
    value_2 = [value_new5;value_new6;value_new7;value_new8]+value1(k);
    value2(:,k) = value_2;
end

[State, State_change1, State_change0] = dictionary();
data_dem = zeros(4,2);

for i = 1:4
    [ row(i) , low(i)] = find(value2 == max(value2(i,:)));
    situation(i) = low(i)+4*(i-1);
    value(i) = value2(row(i),low(i));
    value = value';
    data_dem(i,:) = State(situation(i),:);
    data_dem_1(i,:) = State_change1(situation(i),:);
    data_dem_0(i,:) = State_change0(situation(i),:);
end

value_old = value;
data_dem = [data_dem,zeros(4,data_number-2)]; %前两个时刻的译码值
data_dem_1 = [data_dem_1,zeros(4,data_number-2)]; %多进制分解
data_dem_0 = [data_dem_0,zeros(4,data_number-2)];

for n_dem=3:data_number-D_k(1)
    if mod(n_dem-1,2)+1==1
        for ii=1:2^(L-1)^2*3
            kk=delay+((n_dem-1)*sample_number+1:(n_dem+D_k_new(ii)-1)*sample_number);
            r_k(ii,1)=sum(s_r(kk).*g_kn{ii,1});
        end
    else
        for ii=1:2^(L-1)^2*3
            kk=((n_dem-1)*sample_number+1:(n_dem+D_k_new(ii)-1)*sample_number);
            r_k(ii,1)=sum(s_r(kk).*g_kn{ii,2});
        end
    end
    
    for i = 1:4
        %传3
        sum_index1=sum(data_dem_1(i,:).*h_dem)+1*h(mod(n_dem-1,2)+1);%L = 1
        sum_index0=sum(data_dem_0(i,:).*h_dem)+1*h(mod(n_dem-1,2)+1);%L = 0
        %暴力整伪符号
        b1_0_n = exp(1i*pi*2*(sum_index1));
        b1_0_n_1 = exp(1i*pi*2*(sum_index1-1*h(mod(n_dem-1,2)+1)));
        b1_0_n_2 = exp(1i*pi*2*(sum_index1-1*h(mod(n_dem-1,2)+1)-data_dem_1(i,n_dem-1)*h(mod(n_dem,2)+1)));
        b1_0_n_3 = exp(1i*pi*2*(sum_index1-1*h(mod(n_dem-1,2)+1)-data_dem_1(i,n_dem-1)*h(mod(n_dem,2)+1)-data_dem_1(i,n_dem-2)*h(mod(n_dem-1,2)+1)));
        b1_1_n = exp(1i*pi*2*(sum_index1-data_dem_1(i,n_dem-1)*h(mod(n_dem,2)+1)));
        b1_1_n_1 = exp(1i*pi*2*(sum_index1-1*h(mod(n_dem-1,2)+1)-data_dem_1(i,n_dem-2)*h(mod(n_dem-1,2)+1)));
        b1_2_n = exp(1i*pi*2*(sum_index1-data_dem_1(i,n_dem-2)*h(mod(n_dem-1,2)+1)));
        b1_3_n = exp(1i*pi*2*(sum_index1-data_dem_1(i,n_dem-1)*h(mod(n_dem,2)+1)-data_dem_1(i,n_dem-2)*h(mod(n_dem-1,2)+1)));
        b0_0_n = exp(1i*pi*(sum_index0));
        b0_0_n_1 = exp(1i*pi*(sum_index0-1*h(mod(n_dem-1,2)+1)));
        b0_0_n_2 = exp(1i*pi*(sum_index0-1*h(mod(n_dem-1,2)+1)-data_dem_0(i,n_dem-1)*h(mod(n_dem,2)+1)));
        b0_0_n_3 = exp(1i*pi*(sum_index0-1*h(mod(n_dem-1,2)+1)-data_dem_0(i,n_dem-1)*h(mod(n_dem,2)+1)-data_dem_0(i,n_dem-2)*h(mod(n_dem-1,2)+1)));
        b0_1_n = exp(1i*pi*(sum_index0-data_dem_0(i,n_dem-1)*h(mod(n_dem,2)+1)));
        b0_1_n_1 = exp(1i*pi*(sum_index0-1*h(mod(n_dem-1,2)+1)-data_dem_0(i,n_dem-2)*h(mod(n_dem-1,2)+1)));
        b0_2_n = exp(1i*pi*(sum_index0-data_dem_0(i,n_dem-2)*h(mod(n_dem-1,2)+1)));
        b0_3_n = exp(1i*pi*(sum_index0-data_dem_0(i,n_dem-1)*h(mod(n_dem,2)+1)-data_dem_0(i,n_dem-2)*h(mod(n_dem-1,2)+1)));
        
        [a_k] = calculate( b1_0_n,b1_0_n_1,b1_0_n_2,b1_0_n_3,b1_1_n,b1_1_n_1,b1_2_n,b1_3_n,b0_0_n,b0_0_n_1,b0_0_n_2,b0_0_n_3,b0_1_n,b0_1_n_1,b0_2_n,b0_3_n);
        
        value_new1(i,1) = real(conj(a_k)*r_k);  %传3路径值
    end
    value_new1 = value_new1 + value_old;    %加上之前的路径值
    
    
    for i = 1:4
        %传1
        sum_index1=sum(data_dem_1(i,:).*h_dem)+1*h(mod(n_dem-1,2)+1);%L = 1
        sum_index0=sum(data_dem_0(i,:).*h_dem)-1*h(mod(n_dem-1,2)+1);%L = 0
        %暴力整伪符号
        b1_0_n = exp(1i*pi*2*(sum_index1));
        b1_0_n_1 = exp(1i*pi*2*(sum_index1-1*h(mod(n_dem-1,2)+1)));
        b1_0_n_2 = exp(1i*pi*2*(sum_index1-1*h(mod(n_dem-1,2)+1)-data_dem_1(i,n_dem-1)*h(mod(n_dem,2)+1)));
        b1_0_n_3 = exp(1i*pi*2*(sum_index1-1*h(mod(n_dem-1,2)+1)-data_dem_1(i,n_dem-1)*h(mod(n_dem,2)+1)-data_dem_1(i,n_dem-2)*h(mod(n_dem-1,2)+1)));
        b1_1_n = exp(1i*pi*2*(sum_index1-data_dem_1(i,n_dem-1)*h(mod(n_dem,2)+1)));
        b1_1_n_1 = exp(1i*pi*2*(sum_index1-1*h(mod(n_dem-1,2)+1)-data_dem_1(i,n_dem-2)*h(mod(n_dem-1,2)+1)));
        b1_2_n = exp(1i*pi*2*(sum_index1-data_dem_1(i,n_dem-2)*h(mod(n_dem-1,2)+1)));
        b1_3_n = exp(1i*pi*2*(sum_index1-data_dem_1(i,n_dem-1)*h(mod(n_dem,2)+1)-data_dem_1(i,n_dem-2)*h(mod(n_dem-1,2)+1)));
        b0_0_n = exp(1i*pi*(sum_index0));
        b0_0_n_1 = exp(1i*pi*(sum_index0+1*h(mod(n_dem-1,2)+1)));
        b0_0_n_2 = exp(1i*pi*(sum_index0+1*h(mod(n_dem-1,2)+1)-data_dem_0(i,n_dem-1)*h(mod(n_dem,2)+1)));
        b0_0_n_3 = exp(1i*pi*(sum_index0+1*h(mod(n_dem-1,2)+1)-data_dem_0(i,n_dem-1)*h(mod(n_dem,2)+1)-data_dem_0(i,n_dem-2)*h(mod(n_dem-1,2)+1)));
        b0_1_n = exp(1i*pi*(sum_index0-data_dem_0(i,n_dem-1)*h(mod(n_dem,2)+1)));
        b0_1_n_1 = exp(1i*pi*(sum_index0+1*h(mod(n_dem-1,2)+1)-data_dem_0(i,n_dem-2)*h(mod(n_dem-1,2)+1)));
        b0_2_n = exp(1i*pi*(sum_index0-data_dem_0(i,n_dem-2)*h(mod(n_dem-1,2)+1)));
        b0_3_n = exp(1i*pi*(sum_index0-data_dem_0(i,n_dem-1)*h(mod(n_dem,2)+1)-data_dem_0(i,n_dem-2)*h(mod(n_dem-1,2)+1)));
        
        [a_k] = calculate( b1_0_n,b1_0_n_1,b1_0_n_2,b1_0_n_3,b1_1_n,b1_1_n_1,b1_2_n,b1_3_n,b0_0_n,b0_0_n_1,b0_0_n_2,b0_0_n_3,b0_1_n,b0_1_n_1,b0_2_n,b0_3_n);
        value_new2(i,1) = real(conj(a_k)*r_k);  %传1路径值
    end
    value_new2 = value_new2 + value_old;
    
    for i = 1:4
        sum_index1=sum(data_dem_1(i,:).*h_dem)-1*h(mod(n_dem-1,2)+1);%L = 1
        sum_index0=sum(data_dem_0(i,:).*h_dem)+1*h(mod(n_dem-1,2)+1);%L = 0
        %暴力整伪符号
        b1_0_n = exp(1i*pi*2*(sum_index1));
        b1_0_n_1 = exp(1i*pi*2*(sum_index1+1*h(mod(n_dem-1,2)+1)));
        b1_0_n_2 = exp(1i*pi*2*(sum_index1+1*h(mod(n_dem-1,2)+1)-data_dem_1(i,n_dem-1)*h(mod(n_dem,2)+1)));
        b1_0_n_3 = exp(1i*pi*2*(sum_index1+1*h(mod(n_dem-1,2)+1)-data_dem_1(i,n_dem-1)*h(mod(n_dem,2)+1)-data_dem_1(i,n_dem-2)*h(mod(n_dem-1,2)+1)));
        b1_1_n = exp(1i*pi*2*(sum_index1-data_dem_1(i,n_dem-1)*h(mod(n_dem,2)+1)));
        b1_1_n_1 = exp(1i*pi*2*(sum_index1+1*h(mod(n_dem-1,2)+1)-data_dem_1(i,n_dem-2)*h(mod(n_dem-1,2)+1)));
        b1_2_n = exp(1i*pi*2*(sum_index1-data_dem_1(i,n_dem-2)*h(mod(n_dem-1,2)+1)));
        b1_3_n = exp(1i*pi*2*(sum_index1-data_dem_1(i,n_dem-1)*h(mod(n_dem,2)+1)-data_dem_1(i,n_dem-2)*h(mod(n_dem-1,2)+1)));
        b0_0_n = exp(1i*pi*(sum_index0));
        b0_0_n_1 = exp(1i*pi*(sum_index0-1*h(mod(n_dem-1,2)+1)));
        b0_0_n_2 = exp(1i*pi*(sum_index0-1*h(mod(n_dem-1,2)+1)-data_dem_0(i,n_dem-1)*h(mod(n_dem,2)+1)));
        b0_0_n_3 = exp(1i*pi*(sum_index0-1*h(mod(n_dem-1,2)+1)-data_dem_0(i,n_dem-1)*h(mod(n_dem,2)+1)-data_dem_0(i,n_dem-2)*h(mod(n_dem-1,2)+1)));
        b0_1_n = exp(1i*pi*(sum_index0-data_dem_0(i,n_dem-1)*h(mod(n_dem,2)+1)));
        b0_1_n_1 = exp(1i*pi*(sum_index0-1*h(mod(n_dem-1,2)+1)-data_dem_0(i,n_dem-2)*h(mod(n_dem-1,2)+1)));
        b0_2_n = exp(1i*pi*(sum_index0-data_dem_0(i,n_dem-2)*h(mod(n_dem-1,2)+1)));
        b0_3_n = exp(1i*pi*(sum_index0-data_dem_0(i,n_dem-1)*h(mod(n_dem,2)+1)-data_dem_0(i,n_dem-2)*h(mod(n_dem-1,2)+1)));
        
        [a_k] = calculate( b1_0_n,b1_0_n_1,b1_0_n_2,b1_0_n_3,b1_1_n,b1_1_n_1,b1_2_n,b1_3_n,b0_0_n,b0_0_n_1,b0_0_n_2,b0_0_n_3,b0_1_n,b0_1_n_1,b0_2_n,b0_3_n);
        value_new3(i,1) = real(conj(a_k)*r_k);  %传-1路径值
    end
    value_new3 = value_new3 + value_old;
    
    %传-3
    for i = 1:4
        sum_index1=sum(data_dem_1(i,:).*h_dem)-1*h(mod(n_dem-1,2)+1);%L = 1
        sum_index0=sum(data_dem_0(i,:).*h_dem)-1*h(mod(n_dem-1,2)+1);%L = 0
        %暴力整伪符号
        b1_0_n = exp(1i*pi*2*(sum_index1));
        b1_0_n_1 = exp(1i*pi*2*(sum_index1+1*h(mod(n_dem-1,2)+1)));
        b1_0_n_2 = exp(1i*pi*2*(sum_index1+1*h(mod(n_dem-1,2)+1)-data_dem_1(i,n_dem-1)*h(mod(n_dem,2)+1)));
        b1_0_n_3 = exp(1i*pi*2*(sum_index1+1*h(mod(n_dem-1,2)+1)-data_dem_1(i,n_dem-1)*h(mod(n_dem,2)+1)-data_dem_1(i,n_dem-2)*h(mod(n_dem-1,2)+1)));
        b1_1_n = exp(1i*pi*2*(sum_index1-data_dem_1(i,n_dem-1)*h(mod(n_dem,2)+1)));
        b1_1_n_1 = exp(1i*pi*2*(sum_index1+1*h(mod(n_dem-1,2)+1)-data_dem_1(i,n_dem-2)*h(mod(n_dem-1,2)+1)));
        b1_2_n = exp(1i*pi*2*(sum_index1-data_dem_1(i,n_dem-2)*h(mod(n_dem-1,2)+1)));
        b1_3_n = exp(1i*pi*2*(sum_index1-data_dem_1(i,n_dem-1)*h(mod(n_dem,2)+1)-data_dem_1(i,n_dem-2)*h(mod(n_dem-1,2)+1)));
        b0_0_n = exp(1i*pi*(sum_index0));
        b0_0_n_1 = exp(1i*pi*(sum_index0+1*h(mod(n_dem-1,2)+1)));
        b0_0_n_2 = exp(1i*pi*(sum_index0+1*h(mod(n_dem-1,2)+1)-data_dem_0(i,n_dem-1)*h(mod(n_dem,2)+1)));
        b0_0_n_3 = exp(1i*pi*(sum_index0+1*h(mod(n_dem-1,2)+1)-data_dem_0(i,n_dem-1)*h(mod(n_dem,2)+1)-data_dem_0(i,n_dem-2)*h(mod(n_dem-1,2)+1)));
        b0_1_n = exp(1i*pi*(sum_index0-data_dem_0(i,n_dem-1)*h(mod(n_dem,2)+1)));
        b0_1_n_1 = exp(1i*pi*(sum_index0+1*h(mod(n_dem-1,2)+1)-data_dem_0(i,n_dem-2)*h(mod(n_dem-1,2)+1)));
        b0_2_n = exp(1i*pi*(sum_index0-data_dem_0(i,n_dem-2)*h(mod(n_dem-1,2)+1)));
        b0_3_n = exp(1i*pi*(sum_index0-data_dem_0(i,n_dem-1)*h(mod(n_dem,2)+1)-data_dem_0(i,n_dem-2)*h(mod(n_dem-1,2)+1)));
        
        [a_k] = calculate( b1_0_n,b1_0_n_1,b1_0_n_2,b1_0_n_3,b1_1_n,b1_1_n_1,b1_2_n,b1_3_n,b0_0_n,b0_0_n_1,b0_0_n_2,b0_0_n_3,b0_1_n,b0_1_n_1,b0_2_n,b0_3_n);
        value_new4(i,1)=(real(conj(a_k)*r_k));  %传-3路径值
    end
    value_new4 = value_new4 + value_old;
    value_num = [value_new1 value_new2 value_new3 value_new4];
    
    value_numm{1,n_dem} = value_num; %每时刻的路径值
    
    
    for i = 1:4  %找最优路径 i=1，2，3，4时分别对应传向下一个时刻的+3，+1，-1，-3
        [ row(i) , low(i)] = find(value_num == max(value_num(:,i)));
        value1(i,1) = value_num(row(i),low(i));
    end
    data_dem_past = data_dem;
    data_dem_past_1 = data_dem_1;
    data_dem_past_0 = data_dem_0;
    for i = 1:4   %路径更新
        data_dem(i,:) = data_dem_past(row(1,i),:);
        data_dem_1(i,:) = data_dem_past_1(row(1,i),:);
        data_dem_0(i,:) = data_dem_past_0(row(1,i),:);
        data_dem(i,n_dem) = State(i,1);
        data_dem_1(i,n_dem) = State_change1(i,1);
        data_dem_0(i,n_dem) = State_change0(i,1);
    end
    value_old = value1;
    
end
[row1,~] = find(max(value_old) == value_old);%找剩余4个路径中的最优路径
data_dem11 = data_dem(row1,:); %最终译码值

Q = find(data_dem11 == data);
acc = length(Q)/(data_number-4)*100; %正确率
% figure(5)
% subplot(211),plot(data(1:data_number-4)),xlabel('原始数据序列');
% subplot(212),plot(data_dem11(1:data_number-4)),xlabel('解调后的数据序列');

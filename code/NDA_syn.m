%% --- 环境参数设计 ARTM Tier2信号
clc,clear,close all;

h = [4,5]/16;            %调制指数
L = 3;                   %关联长度
M = 4;                   %进制数
A = 1;                   %幅度
data_number = 1000;      %产生符号数
sample_number = 16;      %每符号采样点数
K = log2(M);             %每符号承载的二进制比特数
% Rb = 1e6;                %信息速率（比特率） bps 1Mbps - 50Mbps，1bps步进 单位bps 每个Ts间传送K个比特 Rb=KRs
Rb = 1000;
Tb = 1/Rb;               %比特间隔
R = Rb/K;               %符号率（码元速率） 单位时间内传送的符号数目 单位Baud
T = 1/R;               %符号间隔 码元周期
fs = R*sample_number;   %采样率 每秒从连续信号中提取并组成离散信号的采样个数
Ts = 1/fs;               %采样周期
% fc = 1e6;                %载波频率
fc = 1000;
shape = 'cos';           %成型波形（rec:矩形，cos:余弦）
N = T/Ts;                %论文中的N，即sample_number
EbN0 = 35;
SNR = EbN0+10*log10(K)-10*log10(0.5*fs/R);   %信噪比
data_ori = randi([0,M-1],1,data_number);     %产生M进制比特
data = 2*(data_ori+1)-1-M;                   %映射为双极性幅度
% fd = 100e3;              %载波频偏
fd = 100;
%% --- 基带信号
% [s_base,fi] = CPM_mod( h,L,sample_number,data_number,T,A,shape,data,fc); %基带信号，相位

if length(h) == 1
    h = [h h];
end
delta_T=T/sample_number;
tr=delta_T:delta_T:L*T;
wave_cos=(1-cos(2*pi*tr/L/T))/2/L/T;        %余弦脉冲形状
wave_rec=repmat(1/2/L/T,1,L*sample_number);  %矩形脉冲形状
temp_wave=zeros(1,data_number*sample_number);
for h_n=1:data_number-1
    hn=h(mod(h_n+1,2)+1);
    temp_wave((h_n-1)*sample_number+1)=hn*data(1,h_n);
end
if findstr(shape, 'cos')
    wave_ori=conv(temp_wave,wave_cos);%成形为余弦
else
    wave_ori=conv(temp_wave,wave_rec); %成形为理想矩形波
end
wave=wave_ori(1:end-L*sample_number+1);
% figure                           %画成形波形
% subplot(211)
% plot(t,wave);axis([0 0.01 -30 30]);title('成型波形');
% xlabel('时间(s)');ylabel('幅度');
temp_q=zeros(1,sample_number*data_number);
temp_q(1)=wave(1)*delta_T;                   %对成型波形做面积积分
for i=2:sample_number*data_number
    temp_q(i)=temp_q(i-1)+wave(i)*delta_T;
    q=temp_q;
end
fi=2*pi*q;               %时变相位，双h已在卷积时带入
s_base=A*exp(1j*fi);     %等效基带信号
% subplot(212)           %画积分波形
% plot(t,q);axis([0 0.2 -3 3]);title('积分波形');xlabel('时间(s)');ylabel('幅度');

% figure
% subplot(211)
% plot(real(s_base)),axis([0 data_number -1.2 1.2]);
% xlabel('信号样本点'),ylabel('幅度/V')
% title('ARTM Tier2信号同相基带波形')
% subplot(212)
% plot(imag(s_base)),axis([0 data_number -1.2 1.2]);
% xlabel('信号样本点'),ylabel('幅度/V')
% title('ARTM Tier2信号正交基带波形')
t=0:1/fs:(length(s_base)-1)/fs;

%% --- 上变频 载波频偏fd
signal_tmp = s_base.*exp(1j*2*pi*fc*t);   %搬移到载波上
signal_fd = signal_tmp.*exp(1j*2*pi*fd*t);    %载波频偏fd
signal = real(signal_fd);
% figure
% plot(imag(fft(signal)))
% figure
% plot(signal),axis([0 100 -1.2 1.2]);
% xlabel('信号样本点'),ylabel('幅度/V')
% title('上变频后ARTM Tier2信号波形')

%% --- AWGN信道
signal_noise = awgn(signal,SNR);
% figure
% plot(signal_noise),axis([0 200 -3 3]);
% xlabel('信号样本点'),ylabel('幅度/V')
% title('加入噪声后信号波形')

%% --- 下变频(通过低通滤波)
[fi_i,fi_q,b] = BaseRecovery(fc,fs,t,signal_noise);   %滤波器
sr = fi_i+1j*fi_q;
% figure
% plot(imag(fft((sr))))
% tst = dfilt.delay(400);
% stt = filter(tst,sr);
% figure 
% plot(real(stt)),axis([0 data_number -1.2 1.2]);
% figure
% subplot(211)
% plot(real(sr)),axis([0 data_number -1.2 1.2]);
% subplot(212)
% plot(imag(sr)),axis([0 data_number -1.2 1.2]);

%% --- 载波频率同步 按照钟声论文的公式

%锁相环参数
C1 = 2^-1;
C2 = 2^-9;
f_e = zeros(1,length(ef));
f_e(1) = 0;
for i = 2:1:length(ef)
    f_e(i) = f_e(i-1)+C1*(ef(i)-ef(i-1))+C2*ef(i);
end
% fdd = zeros(1,length(ef));
% for i = 1:length(ef)
%     fdd(i) = fd;
% end
figure,plot(f_e);    
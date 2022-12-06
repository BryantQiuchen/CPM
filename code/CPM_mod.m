function [s_base,fi] = CPM_mod(h,L,sample_number,data_number,T,A,shape,data,fc)
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
s_base=A*exp(1j*fi);                      %等效基带信号
% subplot(212)                                  %画积分波形
% plot(t,q);axis([0 0.2 -3 3]);title('积分波形');xlabel('时间(s)');ylabel('幅度');
end


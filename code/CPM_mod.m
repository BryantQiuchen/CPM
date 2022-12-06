function [s_base,fi] = CPM_mod(h,L,sample_number,data_number,T,A,shape,data,fc)
if length(h) == 1
    h = [h h];
end
delta_T=T/sample_number;
tr=delta_T:delta_T:L*T;
wave_cos=(1-cos(2*pi*tr/L/T))/2/L/T;        %����������״
wave_rec=repmat(1/2/L/T,1,L*sample_number);  %����������״
temp_wave=zeros(1,data_number*sample_number);
for h_n=1:data_number-1
    hn=h(mod(h_n+1,2)+1);
    temp_wave((h_n-1)*sample_number+1)=hn*data(1,h_n);
end
if findstr(shape, 'cos')
    wave_ori=conv(temp_wave,wave_cos);%����Ϊ����
else
    wave_ori=conv(temp_wave,wave_rec); %����Ϊ������β�
end
wave=wave_ori(1:end-L*sample_number+1);
% figure                           %�����β���
% subplot(211)
% plot(t,wave);axis([0 0.01 -30 30]);title('���Ͳ���');
% xlabel('ʱ��(s)');ylabel('����');
temp_q=zeros(1,sample_number*data_number);
temp_q(1)=wave(1)*delta_T;                   %�Գ��Ͳ������������
for i=2:sample_number*data_number
    temp_q(i)=temp_q(i-1)+wave(i)*delta_T;
    q=temp_q;
end
fi=2*pi*q;               %ʱ����λ��˫h���ھ��ʱ����
s_base=A*exp(1j*fi);                      %��Ч�����ź�
% subplot(212)                                  %�����ֲ���
% plot(t,q);axis([0 0.2 -3 3]);title('���ֲ���');xlabel('ʱ��(s)');ylabel('����');
end


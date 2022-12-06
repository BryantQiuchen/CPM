function [fi_i,fi_q,b] = BaseRecovery(fc,fs,t,sr)
    sr_i = sr.*cos(2*pi*fc*t);
    sr_q = -sr.*sin(2*pi*fc*t);
    b = fir1(60,fc/2/(fs/2));        %设置低通滤波器
    fi_i = 2*filter(b,1,sr_i);       %通过低通滤波器，通过这个滤波器有两个符号延迟
    fi_q = 2*filter(b,1,sr_q);                      
end 








clc,clear;
set(0,'defaultfigurecolor','w');                 %设置图窗背景为白色
EbN0 = 1:1:10;
MLSDAccRate =[58.0315 70.5641 78.5654 90.5654 92.4641 92.5 94.6565 99.01 99.9154 99.9924];
PAMzuiyouAccRate = [57.0835 68.43117 75.52118 89.8955 87.44966 92.5 93.95832 98.54133 99.9 99.9907];
PAMciyouAccRate = [45.3132 70.43117 70.52118 80.8955 82.44966 88.5464 90.4563 93.5461 97.7464 99.8231];
PAMzuiyouErrRate = 1-PAMzuiyouAccRate/100;
MLSDErrRate = 1-MLSDAccRate/100;
PAMciyouErrRate = 1-PAMciyouAccRate/100;
figure
semilogy(EbN0,MLSDErrRate,'bo-',EbN0,PAMzuiyouErrRate,'rd-',EbN0,PAMciyouErrRate,'m>-');
xlabel('Eb/N0(dB)'),ylabel('BER'),legend('MLSD','PAM最优解调','PAM低复杂度解调');
title('误码率'),grid on;
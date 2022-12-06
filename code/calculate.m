%计算伪符号
function [a_k] = calculate( b1_0_n,b1_0_n_1,b1_0_n_2,b1_0_n_3,b1_1_n,b1_1_n_1,b1_2_n,b1_3_n,b0_0_n,b0_0_n_1,b0_0_n_2,b0_0_n_3,b0_1_n,b0_1_n_1,b0_2_n,b0_3_n)
    a_k(1,1) = b0_0_n*b1_0_n;
    a_k(1,2) = b0_0_n_1*b1_0_n;
    a_k(1,3) = b0_0_n*b1_0_n_1;
    a_k(1,4) = b0_0_n_2*b1_0_n;
    a_k(1,5) = b0_0_n*b1_0_n_2;
    a_k(1,6) = b0_0_n_3*b1_0_n;
    a_k(1,7) = b0_0_n*b1_0_n_3;
    a_k(1,8) = b0_1_n*b1_0_n;
    a_k(1,9) = b0_1_n_1*b1_0_n;
    a_k(1,10) = b0_1_n*b1_0_n_1;
    a_k(1,11) = b0_1_n*b1_0_n_2;
    a_k(1,12) = b0_1_n*b1_0_n_3;
    a_k(1,13) = b0_2_n*b1_0_n;
    a_k(1,14) = b0_2_n*b1_0_n_1;
    a_k(1,15) = b0_2_n*b1_0_n_2;
    a_k(1,16) = b0_2_n*b1_0_n_3;
    a_k(1,17) = b0_3_n*b1_0_n;
    a_k(1,18) = b0_3_n*b1_0_n_1;
    a_k(1,19) = b0_3_n*b1_0_n_2;
    a_k(1,20) = b0_3_n*b1_0_n_3;
    a_k(1,21) = b0_0_n*b1_1_n;
    a_k(1,22) = b0_0_n_1*b1_1_n;
    a_k(1,23) = b0_0_n*b1_1_n_1;
    a_k(1,24) = b0_0_n_2*b1_1_n;
    a_k(1,25) = b0_0_n_3*b1_1_n;
    a_k(1,26) = b0_1_n*b1_1_n;
    a_k(1,27) = b0_1_n_1*b1_1_n;
    a_k(1,28) = b0_1_n*b1_1_n_1;
    a_k(1,29) = b0_2_n*b1_1_n;
    a_k(1,30) = b0_2_n*b1_1_n_1;
    a_k(1,31) = b0_3_n*b1_1_n;
    a_k(1,32) = b0_3_n*b1_1_n_1;
    a_k(1,33) = b0_0_n*b1_2_n;
    a_k(1,34) = b0_0_n_1*b1_2_n;
    a_k(1,35) = b0_0_n_2*b1_2_n;
    a_k(1,36) = b0_0_n_3*b1_2_n;
    a_k(1,37) = b0_1_n*b1_2_n;
    a_k(1,38) = b0_1_n_1*b1_2_n;
    a_k(1,39) = b0_2_n*b1_2_n;
    a_k(1,40) = b0_3_n*b1_2_n;
    a_k(1,41) = b0_0_n*b1_3_n;
    a_k(1,42) = b0_0_n_1*b1_3_n;
    a_k(1,43) = b0_0_n_2*b1_3_n;
    a_k(1,44) = b0_0_n_3*b1_3_n;
    a_k(1,45) = b0_1_n*b1_3_n;
    a_k(1,46) = b0_1_n_1*b1_3_n;
    a_k(1,47) = b0_2_n*b1_3_n;
    a_k(1,48) = b0_3_n*b1_3_n;  
%----子函数GetBetakiDk----%
%功能：得到由L产生的各个参数
%输出：二进制分解系数Beta_k_i，波形平移需要的参数矩阵Move_T，k行的C_k的长度矩阵D_k
%输入：关联长度L
function [Beta_k_i,D_k]=GetBetakiDk(L)
Q=2^(L-1);
beta_k_i=zeros(Q,L-1);
D_k=[];
if L~=1
    switch L
        case 2
            beta_k_i=[0;1];
        case 3
            beta_k_i=[0 0;1 0;0 1;1 1];
        case 4
            beta_k_i=[0 0 0;1 0 0;0 1 0;1 1 0;0 0 1;1 0 1;0 1 1;1 1 1];
        case 5
            beta_k_i=[0 0 0 0;1 0 0 0;0 1 0 0;1 1 0 0;0 0 1 0;1 0 1 0;0 1 1 0;1 1 1 0;0 0 0 1;1 0 0 1;0 1 0 1;1 1 0 1;0 0 1 1;1 0 1 1;0 1 1 1;1 1 1 1];
    end
    Beta_k_i=cat(2,zeros(Q,1),beta_k_i);                % i从0到L-1列的beta_k_i
    for i=1:Q
        for j=1:L-1
          temp_D_k(i,j)=L*(2-beta_k_i(i,j))-j;
        end
        D_k(i)=min(temp_D_k(i,:),[],2);            % k行的C_k的长度矩阵，每个长度对应于Move_T的每一行的i个U平移相乘的结果
    end  
end
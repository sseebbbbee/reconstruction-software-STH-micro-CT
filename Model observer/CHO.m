function [ t_s, t_n] = CHO( g_s, g_n, U )
% g_s is the matrix where the columns are the signal image vector
% g_n is the matrix where the columns are the noise image vector
% U is the channel matrix
% t_s is a vector including the test statistics for the signal images
% t_n is a vector including the test statistics for the noise images
v_s = U'*g_s;
v_n = U'*g_n;
K_s=cov(v_s');
K_n=cov(v_n');
K_g=0.5*(K_s+K_n);
g_s_mean=mean(g_s');
g_n_mean=mean(g_n');
delta_v=U'*(g_s_mean'-g_n_mean');
w_CHO=K_g^(-1)*delta_v;
t_s=w_CHO'*v_s;
t_n=w_CHO'*v_n;

end


function D = create_k_mer_D_matrix_noisydup_j_v(k,alphasize,ell,j,v)
M = alphasize^k;
D = zeros(1,M*(M+1)/2);
for m = 1:alphasize-1
    d = create_k_mer_D_matrix_noisydup_j_v_m(k,alphasize,ell,j,v,m);
    D = D+d;
end
D = D/(alphasize-1);
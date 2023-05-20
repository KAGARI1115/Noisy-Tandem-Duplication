function D = create_k_mer_D_matrix_noisydup(k,alphasize,ell)
M = alphasize^k;
N = alphasize^(2*k-2);
D = zeros(M*(M+1)/2,N);
for j = 1:ell
    d = create_k_mer_D_matrix_noisydup_j(k,alphasize,ell,j);
    D = D+d;
end
D = D/ell;
end
    
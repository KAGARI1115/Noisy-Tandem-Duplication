function D= create_k_mer_D_matrix(q,k,alphasize,ell)
Dflip = create_k_mer_D_matrix_noisydup(k,alphasize,ell);
Ddup = create_k_mer_D_matrix_dup(k,alphasize,ell);
D = q(1)*Ddup + q(2)*Dflip;
end
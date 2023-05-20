% Gives the characteristic matrix A for noisy duplication of length ell
% with distribution q, when k>2ell.
function A = create_k_mer_A_matrix_noisydup_largek(k,alph,ell,q)
M = alph^k;
A= zeros(M,M);
for d = 1:ell
    Ad = create_k_mer_A_matrix_noisydup_large_k_new(k,alph,ell,d);
    A = Ad*q(d);
end
end
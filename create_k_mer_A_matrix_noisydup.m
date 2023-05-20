function A = create_k_mer_A_matrix_noisydup(k,alph,ell,q)
if k > 2*ell
    A = create_k_mer_A_matrix_noisydup_largek(k,alph,ell,q);
else
    A = create_k_mer_A_matrix_noisydup_smallk(k,alph,ell,q);
end
end

function A = create_k_mer_A_matrix_noisydup_largek(k,alph,ell,q)
qd = zeros(1,ell+1);
qd(end)=1;
A = q(1)*create_k_mer_A_matrix(qd,k,alph);
for d = 2:ell+1
    Ad = create_k_mer_A_matrix_noisydup_large_k_new(k,alph,ell,d-1);
    A = A+ Ad*q(d);
end
end

function A = create_k_mer_A_matrix_noisydup_smallk(k,alph,ell,q)
qd = zeros(1,ell+1);
qd(end)=1;
A = q(1)*create_k_mer_A_matrix(qd,k,alph);
for d = 2:ell+1
    Ad = create_k_mer_A_matrix_noisydup_small_k_new(k,alph,ell,d-1);
    A = A+Ad*q(d);
end
end
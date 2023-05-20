function P = create_k_mer_D_matrix_noisydup_j_v_m(k,alphasize,ell,j,v,m)
M = alphasize^k;
O = zeros(M,M);
n = 1;
for i = 1:M
    for h = i:M
        O(i,h)=n;
        O(h,i)=n;
        n = n+1;
    end
end
K = diff_k_mer_noisedup(v,k,alphasize,ell,j,m);
P = zeros(1,M*(M+1)/2);
for u = 1:M
    for w = u:M
        n = O(u,w);
        P(n) = K(u)*K(w);
    end
end
end
        
        


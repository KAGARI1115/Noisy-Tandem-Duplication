%%how is the expected difference of second moments depend on the
%%(2k+j-2-ell)-mers after a noisy duplication of length ell with noisy at
%%position j
%%we are expected to see a matrix of size
%%(alphasize)^{2k+j-ell-2}->(alphasize)^{2k-2}
function Dj = create_k_mer_D_matrix_noisydup_j(k,alphasize,ell,j)
M = alphasize^k;
N = alphasize^(2*k+j-2-ell);
D = zeros(M*(M+1)/2,N);
%
O = zeros(M,M);
n = 1;
for i = 1:M
    for h = i:M
        O(i,h)=n;
        O(h,i)=n;
        n = n+1;
    end
end
for vnum = 1:N %%going over all (2k-ell-2)-mers
    v = num2seq(vnum,alphasize,2*k+j-ell-2);
    P = create_k_mer_D_matrix_noisydup_j_v(k,alphasize,ell,j,v)   
    for i = 1:M
        for h = i:M
            n = O(i,h);
            D(n,vnum) = P(n); 
        end
    end
end
%%Then we expand D to D1, D1 is the same thing just of the larger size
Dprime = D.';
expand = alphasize^(ell-j);
Dj = [];
for s = 1:N
    for m = 1:expand
        Dj = [Dj;Dprime(s,:)];
    end
end
Dj = Dj.';
end
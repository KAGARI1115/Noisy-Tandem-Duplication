%%how is the expected difference of second moments depend on the
%%(2k-2-ell)-mers after a duplication of length ell
%%we are expected to see a matrix of size
%%(alphasize)^{2k-ell-2}->(alphasize)^{2k-2}
function D1 = create_k_mer_D_matrix_dup(k,alphasize,ell)
M = alphasize^k;
N = alphasize^(2*k-2-ell);
D = zeros(M*(M+1)/2,N);
%
O = zeros(M,M);
n = 1;
for i = 1:M
    for j = i:M
        O(i,j)=n;
        O(j,i)=n;
        n = n+1;
    end
end
%
for i = 1:M
    for j = i:M
        %%starting to compute the O(i,j)-th row of D
        %%let v be the i-th sequence of k-mers, w be the j-th sequence of
        %%k-mers, then if a duplication is determined by v
        n = O(i,j);
        for vnum = 1:N %%going over all (2k-ell-2)-mers
            v = num2seq(vnum,alphasize,2*k-ell-2);
            %%get the difference vector 
            Diff = diff_k_mer_dup(v,k,alphasize,ell);
            %%get the i-th element of Diff, which is exactly the change
            %%happened to v
            dv = Diff(i);
            %%get the j-th element of Diff, which is the change of w
            dw = Diff(j);
            %%so dvdw is the coefficient of x_n^v, which is the the vnum-th
            %%elemeng of that row
            D(n,vnum) = dv*dw;
        end
    end
end
%%Then we expand D to D1, D1 is the same thing just of the larger size
Dprime = D.';
expand = alphasize^ell;
D1 = [];
for j = 1:N
    for m = 1:expand
        D1 = [D1;Dprime(j,:)];
    end
end
D1 = D1.';
end

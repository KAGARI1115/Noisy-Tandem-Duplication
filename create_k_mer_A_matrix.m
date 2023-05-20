% Gives Matrix A assuming a tandem duplication, substitution system for
% k-mers and alphsize.
% q(l+1) is the prob for dup of length l (or sub for l=0)
% EXAMPLE
% syms alpha
% q = [alpha,1-alpha]
% A = create_k_mer_A_matrix(q,2,2)
function A = create_k_mer_A_matrix(q,k,alphsize)
ell_max = length(q)-1;
M = alphsize^k;
A = zeros(M,M);
for ell = 0:ell_max
    A_ell = create_A_ell_matrix(k,alphsize,ell);
    A = A + q(ell+1)*A_ell;
end
end

function A_ell = create_A_ell_matrix(k,alphsize,ell)
M = alphsize^k;
A_ell = zeros(M,M);
if ell == 0 % Hamming Ball of radius 1
    for j = 1 : M
    v = num2seq(j,alphsize,k);
        for pos = 1:k
            for diff = 1:alphsize-1
                u = v;
                u(pos) = mod(v(pos)+diff,alphsize);
                i = seq2num(u,alphsize);
                A_ell(i,j) = A_ell(i,j)+1/(alphsize-1);
            end
        end
    A_ell(j,j) = A_ell(j,j) - k;
    end
    return
end
for j = 1:M
    v = num2seq(j,alphsize,k);
    vdup = [v(1:ell) v];
    for b = 1:ell
        u = vdup(ell-b+1:ell-b+k);
        i = seq2num(u,alphsize);
        A_ell(i,j) = A_ell(i,j)+1;
    end
    for b = ell+1:k-1
        vdup = [v(1:b) v(b+1-ell:end)];
        u = vdup(1:k);
        i = seq2num(u,alphsize);
        A_ell(i,j) = A_ell(i,j)+1;
    end
    A_ell(j,j) = A_ell(j,j) - (k-1);
end

end
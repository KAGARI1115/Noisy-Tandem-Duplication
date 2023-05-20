% After a duplication of length \ell, what happens to \mu of k-mers
% ell: dup length
% k: k-mer
% a: subsequence of length 2k-ell-2 where dup happens
function K = diff_k_mer_dup(a,k,alphasize,ell)
l = length(a);
if l ~=2*k-ell-2
    disp('input base sequence length wrong')
else
    numkmer=alphasize^k;
    K = zeros(1,numkmer);
    for i = 1:l-k+1
        s = a(i:i+k-1);
        n = seq2num(s,alphasize);
        K(n) = K(n)-1;
    end
    for i = 1:k-1
        s = [a(i:k-1),a(k-ell:k-ell-1+i)];
        n = seq2num(s,alphasize);
        K(n) = K(n) +1;
    end
end

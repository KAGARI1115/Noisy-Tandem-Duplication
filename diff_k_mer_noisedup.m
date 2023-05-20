% After a noisy duplication of length \ell with noise at position i and noise is m, what happens to \mu of k-mers
% ell: dup length
% k: k-mer
% a: subsequence of length 2k+i-ell-2 where dup happens
% i: the position of noise
function K = diff_k_mer_noisedup(a,k,alphasize,ell,i,m)
l = length(a);
if l ~=2*k+i-ell-2
    disp('input base sequence length wrong')
else
    numkmer=alphasize^k;
    K = zeros(1,numkmer);
    %%eliminate case
    for j = 1:k+i-ell-1
        s = a(j:j+k-1);
        n = seq2num(s,alphasize);
        K(n) = K(n)-1;
    end
    %%born case 1: two parts
    for j = 1:i-1
        s1 = a(j:k-1);
        s2 = a(k-ell:k-ell+j-1);
        s=[s1,s2];
        n = seq2num(s,alphasize);
        K(n) = K(n)+1;
    end
    %%born case 2: three parts
    for j = i:k-1
        s1 = a(j:k-1);
        s2 = a(k-ell:k-ell+i-2);
        s3 = a(k-ell+i-1:k-ell+j-1);
        b = s3(1);
        b = mod(b+m,alphasize);
        s = [s1,s2,b,s3(2:end)];
        n = seq2num(s,alphasize);
        K(n) = K(n) + 1;
    end
    %%born case 3: two parts
    for j = k-ell:i+k-ell-1
        s1 = a(j:k-ell+i-1);
        s2 = a(k-ell+i:k+j-1);
        b = s1(end);
        b = mod(b+m,alphasize);
        s = [s1(1:end-1),b,s2];
        n = seq2num(s,alphasize);
        K(n) = K(n) + 1;
    end
    
end
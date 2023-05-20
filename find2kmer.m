function m=find2kmer(seq1,seq2,k,alphasize)
M = alphasize^k;
O = zeros(M,M);
n = 1;
for i = 1:M
    for j = i:M
        O(i,j)=n;
        O(j,i)=n;
        n = n+1;
    end
end
i = seq2num(seq1,alphasize)
j = seq2num(seq2,alphasize)
m = O(i,j);
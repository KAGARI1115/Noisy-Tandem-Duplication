%%we only assume fixed dup length and single position noise
%%more general matrices construction is left to further discussion
% q = [qdup, qflip]
% k: k-mer
% ell: duplication length
function G= create_k_mer_G_matrix(q,k,alphasize,ell)
M = alphasize^k;
%first of all, create a matrix store the order of pairs of k-mers, this
%order can be directly used as index in matlab array, it starts with 1
O = zeros(M,M);
n = 1;
for i = 1:M
    for j = i:M
        O(i,j)=n;
        O(j,i)=n;
        n = n+1;
    end
end
%Then analyze them one by one 
%first case: exact duplication 
G_dup = zeros(M*(M+1)/2);
A_dup = create_k_mer_diff_matrix_dup(k,alphasize,ell);
%second case: flip duplication 
G_flip = zeros(M*(M+1)/2);
A_flip = create_k_mer_A_matrix_noisydup(k,alphasize,ell)+ell*eye(M);
for i = 1:M
    for j = i:M 
        % in this sense, i is the index of the first k-mer
        % j is the index of the second k-mer
        % we are going to construct the O(i,j)-th row of G
        n = O(i,j);
        % then we create a matrix B to reassign the coefficients 
        Bi = create_transfermatrix(k,alphasize,i,O);
        Bj = create_transfermatrix(k,alphasize,j,O);
        G_dup(n,:)=A_dup(i,:)*Bj+A_dup(j,:)*Bi;
        G_flip(n,:) = A_flip(i,:)*Bj + A_flip(j,:)*Bi;
    end
end
G = G_dup*q(1) + G_flip*q(2);
end

function A_ell = create_k_mer_diff_matrix_dup(k,alphasize,ell)
M = alphasize^k;
A_ell = zeros(M);
for j = 1:M
    v = num2seq(j,alphasize,k);
    vdup = [v(1:ell) v];
    for b = 1:ell
        u = vdup(ell-b+1:ell-b+k);
        i = seq2num(u,alphasize);
        A_ell(i,j) = A_ell(i,j)+1;
    end
    for b = ell+1:k-1
        vdup = [v(1:b) v(b+1-ell:end)];
        u = vdup(1:k);
        i = seq2num(u,alphasize);
        A_ell(i,j) = A_ell(i,j)+1;
    end
    A_ell(j,j) = A_ell(j,j) - (k-ell-1);
end
end
function B = create_transfermatrix(k,alphasize,i,O)
%let i,j be the index for k-mers w and u
%then B would be the transfer matrix for computing the coefficients come
%from \mu^w E[d\mu^{u}], in other words, i is outside and j is in the
%bracket
M = alphasize^k;
N = M*(M+1)/2;
B = zeros(M,N);
p = O(i,:);
for j = 1:length(p)
    B(j,p(j))=1;
end
end
    

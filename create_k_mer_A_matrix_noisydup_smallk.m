% Gives Matrix A assuming a tandem duplication with uniform noise 
% k-mers and alphsize.
% assume ell+1<= k <2ell 
function A = create_k_mer_A_matrix_noisydup_smallk(k,alphsize,ell)
M = alphsize^k;
A = zeros(M,M);
for j = 1:ell 
    j
    if j <k+1-ell
        A_j = create_k_mer_A_matrix_noise_j_case1(k,alphsize,ell,j);
        A = A + A_j;
    else
        A_j = create_k_mer_A_matrix_noise_j_case2(k,alphsize,ell,j)
        A = A + A_j;
    end
end
A = A/ell;

    
end

function A = create_k_mer_A_matrix_noise_j_case1(k,alphasize,ell,j)
%%In this case we have j<k+1-ell. So y_b^2 will have three parts and y_b^3 only has one part.
if j<1 || j> ell
    disp('wrong noisy position')
else
    M = alphasize^k;
    A = zeros(M,M);
    for vnum = 1:M
        v = num2seq(vnum,alphasize,k);
        %%y_b^1
        for b = 1:j
            B = Hammingball(v,alphasize,b);
            for h = 1:alphasize-1
                u = B(h,:);
                i = seq2num(u,alphasize);
                A(i,vnum) = A(i,vnum) + 1/(alphasize-1);
            end
        end
        %%y_b^2
        %%part 1
        for b = 1:k-ell
            uprime = [v(ell+1-b:ell),v(1:k-b)];
            B = Hammingball(uprime,alphasize,b+j);
            for h = 1:alphasize-1
                u = B(h,:);
                i = seq2num(u,alphasize);
                A(i,vnum) = A(i,vnum) + 1/(alphasize-1);
            end
        end
        %%part2
        for b = k-ell+1:ell-1
            uprimeh = [v(ell+1-b:ell),v(1:ell-b)];
            uprime = [uprimeh,uprimeh];
            uprime = uprime(1:k);
            B =Hammingball(uprime,alphasize,b+j);
            for h = 1:alphasize-1
                u = B(h,:);
                i = seq2num(u,alphasize);
                A(i,vnum) = A(i,vnum) + 1/(alphasize-1);
            end
        end
        %%part3
        for b = ell:k-j
            uprime = [v(1:b),v(1:b)];
            uprime = uprime(1:k);
            B = Hammingball(uprime,alphasize,b+j);
            for h = 1:alphasize-1
                u = B(h,:);
                i = seq2num(u,alphasize);
                A(i,vnum) = A(i,vnum) + 1/(alphasize-1);
            end
        end
        %%y_b^3
        for b = k-j+1:k-1
            uprime = [v(1:b),v(1:b)];
            u = uprime(1:k);
            i = seq2num(u,alphasize);
            A(i,vnum) = A(i,vnum)+1;
        end
    end
    for i = 1:M
        A(i,i) = A(i,i)-(j+k-1);
    end
end
end
        

function A = create_k_mer_A_matrix_noise_j_case2(k,alphasize,ell,j)
%%In this case we have j>=k+1-ell. So y_b^2 will have two parts and y_b^3 will have two parts.
if j<1 || j> ell
    disp('wrong noisy position')
else
    M = alphasize^k;
    A = zeros(M,M);
    for vnum = 1:M
        v = num2seq(vnum,alphasize,k);
        %%y_b^1
        for b = 1:j
            B = Hammingball(v,alphasize,b);
            for h = 1:alphasize-1
                u = B(h,:);
                i = seq2num(u,alphasize);
                A(i,vnum) = A(i,vnum) + 1/(alphasize-1);
            end
        end
        %%y_b^2
        %%part 1
        for b = 1:k-ell
            uprime = [v(ell+1-b:ell),v(1:k-b)];
            B = Hammingball(uprime,alphasize,b+j);
            for h = 1:alphasize-1
                u = B(h,:);
                i = seq2num(u,alphasize);
                A(i,vnum) = A(i,vnum) + 1/(alphasize-1);
            end
        end
        %%part2
        for b = k-ell+1:k-j
            uprimeh = [v(ell+1-b:ell),v(1:ell-b)];
            uprime = [uprimeh,uprimeh];
            uprime = uprime(1:k);
            B =Hammingball(uprime,alphasize,b+j);
            for h = 1:alphasize-1
                u = B(h,:);
                i = seq2num(u,alphasize);
                A(i,vnum) = A(i,vnum) + 1/(alphasize-1);
            end
        end
        %%y_b^3
        %%part1
        for b = k-j+1:ell-1
            uprime = [v(ell+1-b:ell),v(1:ell-b)];
            u = [uprime,uprime];
            u = u(1:k);
            i = seq2num(u,alphasize);
            A(i,vnum) = A(i,vnum)+1;
        end
        %%part2
        for b = ell:k-1
            uprime = [v(1:b),v(1:b)];
            u = uprime(1:k);
            i = seq2num(u,alphasize);
            A(i,vnum) = A(i,vnum)+1;
        end
    end
    for i = 1:M
        A(i,i) = A(i,i)-(j+k-1);
    end
end
end



function B = Hammingball(s,alphasize,i)
B = [];
for h = 1:alphasize-1
    s1 = s;
    a = s(i);
    a = mod(a+h,alphasize);
    s1(i)=a;
    B = [B;s1];
end
end
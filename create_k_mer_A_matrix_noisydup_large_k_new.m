% Outputs characteristic Matrix $A$ for noisy duplication of length $ell$ and substitution number $d$ 
% The d substitutions are uniformly distributed over the template.
% Inputs: $k$ for k-mers, $alph$ for alphabet size, $ell$ for duplication length, $d$ for number of substitutions.
% This function accounts for the case when k>2ell
% The duplication is characterized as a_{i+1}...a_{i+ell} ->
% a_{i+1}...a_{i+ell}a_{i+1}'...a_{i+ell}', where i denoted the random
% position of the duplication.
% The increased k-mers are denoted y_1,...,y_{k+ell-1}, where y_b denotes
% the increased k-mer ending with y_{i+b} (with or without ').
function A = create_k_mer_A_matrix_noisydup_large_k_new(k,alph,ell,d)
M = alph^k;
A = zeros(M,M);

Epos = subpos(ell,d);
szE = size(Epos);
for j = 1:szE(1)
    A_j = create_k_mer_A_matrix_large_k_noise_fix(k,alph,ell,Epos(j,:));
    A = A + A_j;
end
A = A/szE(1);
end

function A = create_k_mer_A_matrix_large_k_noise_fix(k,alph,ell,evec)
%%In this modulo, we find A given the substitution vector, which indicates
%%the positions where substitutions occur.
M = alph^k;
A = zeros(M,M);
for vnum = 1:M
    v = num2seq(vnum,alph,k);
    %Part 1
    for b = 1:ell
        e = evec(1:b);
        epos = find(e);
        nume = length(epos);
        X = Hamball_pos(v(k-b+1-ell:k-ell),alph,epos);
        for j = 1:(alph-1)^nume
            x = X(j,:);
            u = [v(1:k-b),x];
            ui = seq2num(u,alph);
            A(ui,vnum) = A(ui,vnum) + 1/(alph-1)^nume;
        end
    end
    %Part 2
    for b = ell+1:k-(ell+1)
        e = evec;
        epos = find(e);
        nume = length(epos);
        X = Hamball_pos(v(k-b-ell+1:k-b),alph,epos);
        for i = 1:(alph-1)^nume
            x = X(i,:);
            u = [v(1:k-b),x,v(k-b+1:k-ell)];
            ui = seq2num(u,alph);
            A(ui,vnum) = A(ui,vnum) + 1/(alph-1)^(nume);
        end
    end
    %Part 3
    for b = k-ell:k-1
        e1 = evec(1:ell-(k-b));
        e2 = evec(ell-(k-b)+1:ell);
        epos1 = find(e1);
        epos2 = find(e2);
        nume1 = length(epos1);
        nume2 = length(epos2);
        X = Hamball_pos(v(1:ell-(k-b)),alph,epos1);
        Y = Hamball_pos(v(ell-(k-b)+1:ell),alph,epos2);
        for i = 1:(alph-1)^nume1
            for j = 1:(alph-1)^nume2
                x = X(i,:);
                y = Y(j,:);
                u = [v(ell-(k-b)+1:ell),x,y,v(ell+1:b)];
                ui = seq2num(u,alph);
                A(ui,vnum) = A(ui,vnum) + 1/(alph-1)^(nume1+nume2);
            end
        end
    end
    %Part 4
    for b = k:k+ell-1
        e = evec(b+1-k:ell);
        epos = find(e);
        nume = length(epos);
        X = Hamball_pos(v(1:ell+k-b),alph,epos);
        for i = 1:(alph-1)^nume
            x = X(i,:);
            u = [x, v(ell+k-b+1:k)];
            ui = seq2num(u,alph);
            A(ui,vnum) = A(ui,vnum) + 1/(alph-1)^(nume);
        end
    end 
end
for i = 1:M
    A(i,i) = A(i,i)-(ell+k-1);
end
end
% This function produces a matrix of size {ell\choose d} \times ell, where each row contains a binary vector.
% Each binary vector corresponds to an instance of choosing d from ell, the
% 1's indicates the chosen items.
function Epos = subpos(ell,d)
    if d == 1
        Epos = zeros(ell,ell);
        for j = 1:ell
            Epos(j,j) = 1;
        end
    elseif d==0
        Epos = [];
    else
        Epos1 = subpos(ell,d-1);
        Epos = zeros(nchoosek(ell,d),ell);
        sz = size(Epos1);
        k = 1;
    for i = 1:sz(1)
        vec = Epos1(i,:);
        j = length(vec);
        while vec(j) == 0
            vecnew = vec;
            vecnew(j)=1;
            Epos(k,:) = vecnew;
            j = j-1;
            k = k+1;
        end
    end


    end
end
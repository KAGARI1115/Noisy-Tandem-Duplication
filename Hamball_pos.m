function Ham = Hamball_pos(s,alphasize,pos)
if length(pos) == 1
    p = pos(1);
    Ham = zeros(alphasize-1, length(s));
    for j = 1:alphasize-1
        mu = s;
        mu(p) = mod(mu(p)+ j,alphasize);
        Ham(j,:) =mu;
    end
elseif isempty(pos)
    Ham = s;
else
    Ham1 = Hamball_pos(s,alphasize,pos(2:end));
    sz = size(Ham1);
    Ham = zeros(sz(1)*(alphasize-1), length(s));
    k = 1;
    p = pos(1)
    for i = 1:sz(1)
        base = Ham1(i,:);
        for j = 1:alphasize-1
            mu = base;
            mu(p) = mod(mu(p)+ j,alphasize);
            Ham(k,:) = mu;
            k = k+1;
        end
    end
end
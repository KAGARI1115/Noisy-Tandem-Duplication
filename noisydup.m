function s =noisydup(s0,alphasize,ell,d) %Basically it is for binary and dup length 1
l = length(s0);
if l >= ell
    s = [s0,s0];
    i = randi(l);
    tem = s(i:i+ell-1);
    E = subpos(ell,d);
    sz = size(E);
    j = randi([1,sz(1)]);
    sub = E(j,:);
    subp = find(sub);
    for k = 1:length(subp)
        a = randi(alphasize-1);
        pp = subp(k);
        tem(pp) = mod(tem(pp) + a,alphasize);
    end
    s = [s(1:i+ell-1),tem,s(i+ell:end)];
    s = s(i:i+ell+l-1);
end
end

function X = countfre(s,k,alphasize)
X = zeros(alphasize^k,1);
l =length(s);
s = [s,s(1:k-1)];
for i = 1:l
    a = stringtonum(s(i:i+k-1),alphasize);
    X(a+1) = X(a+1)+1;
end
X = X/l;
end



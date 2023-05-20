function a = stringtonum(s,alphasize)
l = length(s);
a = 0;
for i = l:-1:1
    a = a+ s(i)*(alphasize^(l-i));
end


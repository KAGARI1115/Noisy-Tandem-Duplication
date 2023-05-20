function s =flipdup(A,alphasize) %Basically it is for binary and dup length 1
l = length(A);
s = A;
if l >=1
i = randi(l);
a = randi(alphasize-1);
a = mod(s(i) + a,alphasize);
s = [s(1:i), a, s(i+1:end)];
end
end

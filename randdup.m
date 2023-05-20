function s =randdup(A,d)
l = length(A);
if l < d
    s = A;
else
A = [A(end-d+2:end),A];
i = randi(l);
if i+d<=length(A)
    s = [A(1:i+d-1),A(i:i+d-1),A(i+d:end)];
else
    s = [A(1:i+d-1),A(i:i+d-1)];
end
s(1:d-1) = [];
end


    
function B=randdelete(A)
l = length(A);
if l == 0
    B=A;
else
i = randi(l);
B = A;
B(i) = [];
end
end

function B = randsub(A)
B = A;
l = length(A);
if l ~= 0
    
i = randi(l);
a = rand;
if a<=0.333
    b = B(i)-1;
    b = mod(b+1,4);
    B(i) = b+1;
elseif a<=0.667
    b = B(i)-1;
    b = mod(b+2,4);
    B(i) = b+1;
else
    b = B(i)-1;
    b = mod(b+1,4);
    B(i) = b+1;
end
end
end
    
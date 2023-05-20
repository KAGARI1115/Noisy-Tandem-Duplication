%insert into position between and at the end, random
function B = randinsert(A)
l = length(A);
if l <= 1
    B=A;
else
i = randi(l-1);
a = randi(4);
B = [A(1:i),a,A(i+1:end)];
end
end



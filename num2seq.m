% convert a number to the sequence of length k whose lexigraphic rank is
% num. 
% elements of the sequence can be 0, 1, ..., alphsize-1
% the LSB is the largest index in seq
% example:
% num2seq(4,2,5)
% ans =
%      0     0     0     1     1
function seq = num2seq(num,alphsize,k)
num = num - 1;
seq = [];
for i = 1 : k
    seq = [mod(num,alphsize), seq]; %#ok<AGROW>
    num = floor(num/alphsize);
end
end
% give the lexiographic rank of the sequences (rank starts from 1)
% elements of the sequence can be 0, 1, ..., alphsize-1
% the LSB is the largest index in seq
% seq2num([0,0,0,1,1],2)
% ans =
%      4
function num = seq2num(seq,alphsize)
num = 0;
for i = 1 : length(seq)
    num = alphsize * num + seq(i);
end
num = num + 1;
end
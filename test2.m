%test A matrix 
syms d
k = 2;
alph= 3;
ell=1;
q(d)=[1-d,d];
Aflip = create_k_mer_A_matrix_noisydup(k,alph,ell,q(d));
[W,D] = eig(Aflip)
inv(W)*[1;0;0;0;0;0;0;0;0]

%% 
s = [0,1,2,3]
B = [];
for h = 1:alphasize-1
    s1 = s;
    a = s(i);
    a = mod(a+h,alphasize);
    s1(i)=a;
    B = [B;s1];
end
%%
a=0.01
(2-a)/(3*a+2)

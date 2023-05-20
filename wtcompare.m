alph = 3;
k=2;
ell=1;
L = 10;
M = alph^k;
Aflip = create_k_mer_A_matrix_noisydup(k,alph,ell);
Adup = create_k_mer_A_matrix([0,1],k,alph);
s0 = zeros(1,L);
N1 = [];
q = [0.05:0.01:0.2];
for d = 1:length(q)
    delta = q(d);
    A = Adup*(1-delta) + Aflip*delta;
    n = 0;
    x0 = countfre(s0,k,alph);
    mu = x0*L;
    while mu(5)<1
        n = n+1;
        x0 = (eye(M)+A/(L+n))*x0;
        mu = x0*(L+n);
    end
    N1 = [N1,n];
end

%%
s1 = zeros(1,L);
v =5;

Rwt = [];
parfor k = 1:length(q)
    k
   Nd=0;
   for j = 1:n
   Nd = Nd+wtime(s1,q(k),2,v,alph);
   end
   Rwt = [Rwt,Nd/n];
   
end
%%
semilogy(q,rwt,'d-')
hold on 
semilogy(q,N,'o-')
hold on 
semilogy(q,Rwt,'*-')
hold on 
semilogy(q,N1,'v-')
ylabel('$n$','interpreter','latex')
xlabel('$\delta$','interpreter','latex')
legend({'expected waiting time of 12','$\hat n_{12}$','expected waiting time of 11','$\hat n_{11}$'},'interpreter','latex')
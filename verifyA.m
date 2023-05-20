alph=2;
k=4;
ell=3;
q = rand(1,ell+1);
q=q/sum(q);

A = create_k_mer_A_matrix_noisydup(k,alph,ell,q);
s0 = randi([0,alph-1],1,10);
x0 = countfre(s0,k,alph);
X_anl = [];
kmer = [0,1,1,0];
m = seq2num(kmer,alph);
L0= length(s0);
nmut=100;
xfre=x0;
L = L0;
for j = 1:nmut
    L = L+ell;
    xfre = xfre + A*xfre/L;
    X_anl = [X_anl,xfre(m)];
end

sample=2000;
X = zeros(sample,nmut);
for i = 1:sample
    i
    s = s0;
    for j = 1:nmut
        a = mnrnd(1,q);
        d = find(a);
        if d == 1
            s = randdup(s,ell);
        else
            s = noisydup(s,alph,ell,d-1);
        end
        x = countfre(s,k,alph);
        X(i,j) = x(m);
    end
end
E = mean(X);
%figure
%plot(1:dp,VV)
%%
plot(1:nmut,E)
hold on 
plot(1:nmut,X_anl)
%%



alph=3;
k=2;
M = alph^k;
L=10;
dlt=0.1;
q = [1-dlt,dlt];
ell=1;
D = create_k_mer_D_matrix(q,k,alph,ell);
G = create_k_mer_G_matrix(q,k,alph,ell);

var0 = zeros(M*(M+1)/2,1);
var0(1)=1;

BX=[];
Var=[];
kprime = 2*k-2;
Aflip = create_k_mer_A_matrix_noisydup(k,alph,ell);
Adup = create_k_mer_A_matrix([0,1],k,alph);
A = q(1)*Adup + q(2)*Aflip;
x0 = countfre(s0,k,alph);
Bflip = create_k_mer_A_matrix_noisydup(kprime,alph,ell);
Bdup = create_k_mer_A_matrix([0,1],kprime,alph);
B = q(1)*Bdup + q(2)*Bflip;
x1 = countfre(s0,kprime,alph);
for i= 1:nmut
    var0 = var0*(L+i-1)^2/(L+i)^2 + G*(L+i-1)/(L+i)^2 * var0+D/((L+i)^2)*x1;
    x1 = (eye(alph^kprime)+B/(L+i))*x1;
    x0 = (eye(alph^k)+A/(L+i))*x0;
    BX= [BX,x0(m)];
    %var0 = (eye(4)+A/(10+i))*var0;
    Var = [Var,var0];
end
 VVV = zeros(1,nmut);
 P1=zeros(1,nmut);
 P2 = zeros(1,nmut);
 O = zeros(M,M);
n = 1;
for i = 1:M
    for h = i:M
        O(i,h)=n;
        O(h,i)=n;
        n = n+1;
    end
end
n = O(m,m);
P1 = [];
P2 = [];
 for k = 1:nmut
     comvar = (Var(n,k))-BX(k)^2;
     VVV(k) = comvar;
     P1 = [P1,BX(k)+3*sqrt(comvar)];
     P2 = [P2,BX(k)-3*sqrt(comvar)];
 end
 
figure()
plot(1:nmut,VV,'*-')
hold on
plot(1:nmut,VVV,'d-')
legend('Variance of $x_n^{12}$ over 10000 trials','$E[(x_n^{12})^2]-E[x_n^{12}]^2$')
xlabel('$n$','interpreter','latex')
set(legend,'Interpreter','latex')

% plot(1:dp,BX)
% hold on 
% plot(1:dp,E)
% legend('$E[\mu_n^{01}]$','Expected value of $\mu_n^{01}$ over 5000 trials')
% set(legend,'Interpreter','latex')

%%
figure()
plot(1:nmut,P1,'--')
hold on 
plot(1:nmut,BX)
hold on 
plot(1:nmut,P2,'-.')
legend({'$\mathbf{E}[\mu_n^{12}]+3\sigma_n^{12}$','$\mathbf{E}[\mu_n^{12}]$','$\mathbf{E}[\mu_n^{12}]-3\sigma_n^{12}$'},'interpreter','latex')
xlabel('$n$','interpreter','latex')
%%
n=4500;
P =[];
m=1;
p=0;
while p>=0
    %m = BX(n)-gamma*VVV(n);
    %M = [M,m];
    gamma =(BX(n)-m)/VVV(n); 
    p = 1-1/gamma^2;
    P = [P,p];
    m = m+1;
end
plot(1:length(P)-1,P(1:end-1))
xlabel('$m$','interpreter','latex')
ylabel('$1-\frac{1}{{\gamma''}^2}$','interpreter','latex')
title('lower bound of $\mathbf{P}(T_{12}(m)<4500)$ vs $m$','interpreter','latex')
%%
figure()
plot(1:nmut,VV,'*-')
hold on
plot(1:nmut,VVV,'d-')
legend('Variance of $x_n^{00}$ over 5000 trials','$E[(x_n^{00})^2]-E[x_n^{00}]^2$')
xlabel('$n$','interpreter','latex')
set(legend,'Interpreter','latex')
 %%

plot(1:nmut,VVV)
hold on 
plot(1:nmut,BX)
% legend('$E[\mu_n^{01}]$','Expected value of $\mu_n^{01}$ over 5000 trials')
% set(legend,'Interpreter','latex')
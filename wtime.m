function Nd = wtime(s,q,k,v,alphasize)
itr = 200;
Nd=0;
for i = 1:itr
    s1 = s;
    N = 0;
    x = countfre(s1,k,alphasize);
    while x(v)==0
        N = N+1;
        a = rand;
        if a<q
            s1 = flipdup(s1,alphasize);
        else 
            s1 = randdup(s1,1);
        end
        x = countfre(s1,k,alphasize);
    end
    Nd = Nd+N;
end
Nd = Nd/itr;
        
        
    

function raiz = resolver(a,b,c)
%resolve equaçao do 2 grau e retorna valor positivo
    raiz=[];
    delta=b^2-4*a*c;
    raiz1=(-b+sqrt(delta))/(2*a);
    raiz2=(-b-sqrt(delta))/(2*a);
    if raiz1>0
        raiz=raiz1;
    end
    if raiz2>0
        raiz=[raiz, raiz2];
    end
end


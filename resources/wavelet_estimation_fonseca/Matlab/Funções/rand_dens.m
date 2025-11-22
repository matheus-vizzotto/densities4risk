function [ g ] = rand_dens(x,j0,j1,n,wname,cs_phi,cs_psi)
% fornecidos os níveis j0 e j1 da decomposição de ondaletas, os valores dos
% coeficientes e o n da aproximação por Daubechies-Lagarias, essa função
% retorna o valor da densidade resultante avaliada no ponto x.
% Aqui a densidade é chamada aleatória porque supomos que os coeficientes
% foram gerados aleatoriamente.

h = fliplr(wfilters(wname));
N = length(h)/2;
    
gs = 0;
for k = (1-2*N):(2^j0)
    gs = gs + cs_phi(k+2*N)*(2^(j0/2))*phi_n((2^j0)*x-k,n,wname);
end

ki = 1;
for jj = j0:j1
    for k = (-N):(2^jj + N-1)
        gs = gs + cs_psi(ki)*(2^(jj/2))*psi_n((2^jj)*x-k,n,wname);
        ki = ki + 1;
    end
end
g = gs^2;

end


function [ y ] = FuncMatProdYcentEigFunc( fh_Ycent, fh_OrtNormEigFunc, x )
% Essa função retorna uma matriz nxd0 contendo o produto dos funcionais
% observados centrados e das autofunções do operador, ambos avaliados em x,
% que deve ser um escalar.

% tamanho amostral
n = length(fh_Ycent);
% número de autofunções (dimensão do funcional sob H0)
d0 = length(fh_OrtNormEigFunc);
% vetor para armazenar os funcionais observados centrados avaliados em x
vYcent = zeros(n,1);
% vetor para armazenar as autofunções avaliadas em x
vEigFunc = zeros(d0,1);

for ii=1:d0
    vEigFunc(ii) = fh_OrtNormEigFunc{ii}(x);
end

for ii=1:n
    vYcent(ii) = fh_Ycent{ii}(x);
end

% matriz com o produto entre autofunções e funcionais centrados
y = vYcent*vEigFunc';

end


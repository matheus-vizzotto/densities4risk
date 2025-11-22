function [ y ] = FuncMatProdYcent( fh_Y, x)
% A função retorna uma matriz de dimensão nxn com o produto dos funcionais
% observados na célula fh_Y centrados pela média e avaliados em x, que deve
% ser um escalar.

% tamanho da amostra
n = length(fh_Y);
% vetor que armazenará os funcionais observados centrados
vY = zeros(n,1);

for ii=1:n
    vY(ii) = fh_Y{ii}(x);
end
vY = vY - mean(vY);

y = vY*vY';

end


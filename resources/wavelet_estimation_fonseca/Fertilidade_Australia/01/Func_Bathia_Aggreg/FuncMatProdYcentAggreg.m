function [ y ] = FuncMatProdYcentAggreg( fh_Y, w, p, k, x)
% A função retorna uma matriz de dimensão (n-p-delta+1)x(n-p-delta+1) com o
% produto dos funcionais observados na célula fh_Y centrados pela média e 
% multiplicados pelos pesos no vetor w. As funções são avaliadas em x, que deve
% ser um escalar. O lag máximo do processo é dado em p e o lag usado com os
% funcionais ponderados são armazenados em k.

% tamanho da amostra
n = length(fh_Y);
% vetor que armazenará os funcionais observados centrados
vY = zeros(n,1);
% número de pesos no vetor w
delta = length(w);

% calculando os funcionais observados no ponto x
for ii=1:n
    vY(ii) = fh_Y{ii}(x);
end
vY = vY - mean(vY);
% obtendo os funcionais observados ponderados
mY = zeros(delta,delta);
vZ = zeros(n-p-delta+1,1);
for jj = 1:(n-p-delta+1)
    for ii = 0:(delta-1)
        mY(ii+1,:) = vY((jj+k-ii):(jj+k+delta-1-ii))';
    end
    vZ(jj) = w'*mY*w;
end

% matriz com produtos dos funcionais observados ponderados
y = vZ*vZ';

end


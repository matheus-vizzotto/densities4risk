function [ y ] = FuncYcent( fh, ind, n, x)
% função que retorna o funcional ind centrado pela média. fh é uma célula
% que armazena os funcionais, de uma amostra de tamanho n e avaliados em x.

% array que armazenará o funcional centrado
vY = zeros([size(x) n]);

for ii=1:n
    % matrizes dos funcionais em fh avaliados em x
    vY(:,:,ii) = fh{ii}(x);
end

y = fh{ind}(x) -  mean(vY,3);
end


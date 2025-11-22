function [ y ] = FuncLinComb( fh, vec, x)
% a função retorna a combinação linear das funções armazenadas na célula fh
% com os coeficientes do vetor vec. As funções são avaliadas em x.
n = length(vec);
y = 0;
for ii=1:n
    y = y + vec(ii)*fh{ii}(x);
end

% tentativa vetorizada
% n = length(vec);
% vY = zeros([size(x) n]);
% for ii=1:n
%     vY(:,:,ii) = vec(ii)*fh{ii}(x);
% end
% y = sum(vY,3);

end


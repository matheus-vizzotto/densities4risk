function [ y ] = FuncYbar( fh, n, x)
% Função para calcular a função média dos elementos da célula fh.
% n é o tamanho da amostra e x o valor onde as funções serão avaliadas
y=0;
for ii=1:n
    y = y + fh{ii}(x);
end
y = y/n;
end


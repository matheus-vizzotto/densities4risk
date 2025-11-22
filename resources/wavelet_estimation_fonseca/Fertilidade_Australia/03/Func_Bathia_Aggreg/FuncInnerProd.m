function y = FuncInnerProd( func1, func2, xmin, xmax )
% Função para calcular o produto interno entre as funções func1 e func2,
% avaliando a integral numericamento no intervalo de xmin até xmax

func3 = @(x) func1(x).*func2(x);
y = integral(func3,xmin,xmax);

end


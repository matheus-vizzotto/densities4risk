function [ y ] = VecFunc( func, x )
% função usada para vetorizar func, avaliando func em cada elemento da
% matriz x
[n,m] = size(x);
y = zeros(n,1);
for ii=1:n
    for jj=1:m
        y(ii,jj) = func(x(ii,jj));
    end
end

end


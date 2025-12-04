function [ Yort ] = MatGramSchm(Y)
% função para aplicar o método de Gram-Schmidt aos vetores da matriz Y

[n,d] = size(Y);
% matriz para armazenar os vetores ortonormalizados
Yort = zeros(n,d);
Yort(:,1) = Y(:,1);

% etapa de ortogonalização
for ii=2:d
    Yort(:,ii) = Y(:,ii);
    for jj=1:(ii-1)
        Yort(:,ii) = Yort(:,ii) - ((Yort(:,jj)'*Y(:,ii))/(Yort(:,jj)'*Yort(:,jj)))*Yort(:,jj);
    end
end

% etapa de ortonormalização
for ii=1:d
    Yort(:,ii) = Yort(:,ii)/norm(Yort(:,ii));
end

end


function [ vecd0 ] = BootFuncDim( fh_Y, fh_Ycent, EigVecK, mEta, d0, p, Nboot )
% Função para obter uma amostra bootstrap dos autovalores da matriz K* de
% produtos internos. Considerando que se dejesa testar se o (d0+1)-ésimo
% autovalor é zero, serão gerados funcionais bootstrap com dimensão d0,
% para os quais serão calculados o (d0+1)-ésimo autovalor da K*
% correspondente, que serão armazenados em vecd0.

% tamanho da amostra
n = length(fh_Y);
% vetor que armazenará a amostra bootstrap do (d0+1)-ésimo autovalor
vecd0 = zeros(Nboot,1);
% d0 primeiras autofunções do operador
fh_OrtNormEigFunc = FuncGramSchm( fh_Ycent, EigVecK, d0);

for jj=1:Nboot
    % gerando uma amostra bootstrap
    
    % sorteando os índices do resíduo bootstrap reamostrado
    indboot = unidrnd(n,n,1);
    
    % funcional observado centrado
    fh_Ycentboot = cell(n,1);
    for ii=1:n
        fh_Ycentboot{ii} = @(z) FuncYcentboot( fh_Y, fh_OrtNormEigFunc, mEta, indboot, ii, z);
    end

    Kstar = zeros(n-p,n-p);
    for k=1:p
        Gk = zeros(n-p,n-p);
        for tt=1:(n-p)
             for ss=tt:(n-p)
                 innprodGk = InnerProdFunc(fh_Ycentboot{tt+k},fh_Ycentboot{ss+k},0,1);
                 Gk(tt,ss) = innprodGk;
                 Gk(ss,tt) = innprodGk;
             end
        end

        Kstar = Kstar + Gk;
    end
    
    G0 = zeros(n-p,n-p);
    for tt=1:(n-p)
         for ss=tt:(n-p)
             innprodGk = InnerProdFunc(fh_Ycentboot{tt},fh_Ycentboot{ss},0,1);
             G0(tt,ss) = innprodGk;
             G0(ss,tt) = innprodGk;
         end
    end

    Kstar = Kstar*G0/((n-p)^2);
    [~,EigValK] = eig(Kstar);
    vecd0(jj) = EigValK(d0+1,d0+1);
end

end


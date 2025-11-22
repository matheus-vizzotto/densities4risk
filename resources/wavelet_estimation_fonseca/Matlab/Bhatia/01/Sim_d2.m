
% Cenário de simulação usado por Bathia et al (2010) com dimensão 2

n = 50; % tamanho amostral
d = 2; % parâmetro de dimensão
p = 5;

% processos para gerar a dependência temporal da função não observada
% Quando d=2, os coefs são -.65 e .4
model_xi1 = arima('Constant',0,'AR',{-.65},'Variance',1.5);
model_xi2 = arima('Constant',0,'AR',{.4},'Variance',1.5);

xi1 = simulate(model_xi1,n);xi2 = simulate(model_xi2,n);

% funcional observado
fh_Y = cell(n,1);

for ii=1:n
    vZ = normrnd(0,1,10,1);
    fh_Y{ii} = @(x) xi1(ii)*sqrt(2)*cos(pi*x) + xi2(ii)*sqrt(2)*cos(2*pi*x)...
        + sqrt(2)*(vZ(1)*sin(pi*x) + vZ(2)*sin(2*pi*x)/2 ...
        + vZ(3)*sin(3*pi*x)/4 + vZ(4)*sin(4*pi*x)/8 ...
        + vZ(5)*sin(5*pi*x)/16 + vZ(6)*sin(6*pi*x)/32 ...
        + vZ(7)*sin(7*pi*x)/64 + vZ(8)*sin(8*pi*x)/128 ...
        + vZ(9)*sin(9*pi*x)/256 + vZ(10)*sin(10*pi*x))/512;
end

%%%%%%%%%%%%%%%%

% funcional observado centrado
fh_Ycent = cell(n,1);

for ii=1:n
    fh_Ycent{ii} = @(z) fh_Y{ii}(z) - FuncYbar(fh_Y,n,z);
end

Kstar = zeros(n-p,n-p);
for k=1:p
    Gk = zeros(n-p,n-p);
    for tt=1:(n-p)
         for ss=tt:(n-p)
             innprodGk = InnerProdFunc(fh_Ycent{tt+k},fh_Ycent{ss+k},0,1);
             Gk(tt,ss) = innprodGk;
             Gk(ss,tt) = innprodGk;
         end
    end
    
    Kstar = Kstar + Gk;
end

G0 = zeros(n-p,n-p);
for tt=1:(n-p)
     for ss=tt:(n-p)
         innprodGk = InnerProdFunc(fh_Ycent{tt},fh_Ycent{ss},0,1);
         G0(tt,ss) = innprodGk;
         G0(ss,tt) = innprodGk;
     end
end

Kstar = Kstar*G0/((n-p)^2);
[EigVecK,EigValK] = eig(Kstar);

d0 = 2;
mEta = zeros(n,d0);
for ii=1:n
   for jj=1:d0
       mEta(ii,jj) = InnerProdFunc(fh_Ycent{ii},fh_OrtNormEigFunc{jj},0,1);
   end
end


tic;
vecd0 = BootFuncDim( fh_Y, fh_Ycent, EigVecK, 1, p, 2);
toc;



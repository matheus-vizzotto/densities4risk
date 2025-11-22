
% Cenário de simulação usado por Bathia et al (2010) com dimensão 4

n = 20; % tamanho amostral
d = 4; % parâmetro de dimensão
p = 5;

% processos para gerar a dependência temporal da função não observada
% Quando d=4 os coeficientes são -.775,.65,-.525 e .4,
model_xi1 = arima('Constant',0,'AR',{-.775},'Variance',1.5);
model_xi2 = arima('Constant',0,'AR',{.65},'Variance',1.5);
model_xi3 = arima('Constant',0,'AR',{-.525},'Variance',1.5);
model_xi4 = arima('Constant',0,'AR',{.4},'Variance',1.5);

xi1 = simulate(model_xi1,n);xi2 = simulate(model_xi2,n);
xi3 = simulate(model_xi3,n);xi4 = simulate(model_xi4,n);

% funcional observado
fh_Y = cell(n,1);

for ii=1:n
    vZ = normrnd(0,1,10,1);
    fh_Y{ii} = @(x) xi1(ii)*sqrt(2)*cos(pi*x) + xi2(ii)*sqrt(2)*cos(2*pi*x)...
        + xi3(ii)*sqrt(2)*cos(3*pi*x) + xi4(ii)*sqrt(2)*cos(4*pi*x)...
        + sqrt(2)*(vZ(1)*sin(pi*x) + vZ(2)*sin(2*pi*x)/2 ...
        + vZ(3)*sin(3*pi*x)/4 + vZ(4)*sin(4*pi*x)/8 ...
        + vZ(5)*sin(5*pi*x)/16 + vZ(6)*sin(6*pi*x)/32 ...
        + vZ(7)*sin(7*pi*x)/64 + vZ(8)*sin(8*pi*x)/128 ...
        + vZ(9)*sin(9*pi*x)/256 + vZ(10)*sin(10*pi*x))/512;
end

%%%%%%%%%%%%%%%%

% célula para o funcional observado centrado
fh_Ycent = cell(n,1);
for ii=1:n
    fh_Ycent{ii} = @(z) fh_Y{ii}(z) - FuncYbar(fh_Y,n,z);
end

% funcional observado centrado
mat_G = integral(@(x) FuncMatProdYcent( fh_Y, x),0,1,'ArrayValued',true);

Gk = zeros(n-p,n-p,p);
for k=1:p
    Gk(:,:,k) = mat_G((1+k):(n-p+k),(1+k):(n-p+k));
end

Kstar = sum(Gk,3)*mat_G(1:(n-p),1:(n-p))/((n-p)^2);

[EigVecK,EigValK] = eig(Kstar);
diag(EigValK(1:10,1:10))'

% calculando as d0 primeiras autofunções do operador analisado
dmax = 4;
fh_OrtNormEigFunc = FuncGramSchm( fh_Ycent, EigVecK, dmax);
mEta = integral(@(x) FuncMatProdYcentEigFunc( fh_Ycent, fh_OrtNormEigFunc, x ),0,1,'ArrayValued',true);


%tic;
%vecd0 = FuncBootFuncDim(fh_Y, fh_OrtNormEigFunc, mEta, d0, p, 3);
%toc;

rng(2018)
Nboot = 20; alpha = .05; vDimSel = 1;
checkDimSel = 0;
d0 = 1;
mPvalues = ones(1,dmax);
tic;
while d0<=dmax
    fh_OrtNormEigFunc = FuncGramSchm( fh_Ycent, EigVecK, d0);
    d_boot = FuncBootFuncDim(fh_Y, fh_OrtNormEigFunc, mEta(:,1:d0), d0, p, Nboot);
    mPvalues(d0) = sum(d_boot>EigValK(d0+1,d0+1))/Nboot;
    if (mPvalues(d0)>alpha)&&(checkDimSel==0)
        vDimSel = d0;
        checkDimSel = 1;
        %d0 = 7;
    end
    d0 = d0 + 1;
end
toc;
vDimSel
mPvalues


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




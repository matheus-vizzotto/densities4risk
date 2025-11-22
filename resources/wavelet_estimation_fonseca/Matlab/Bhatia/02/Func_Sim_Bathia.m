function [ vDimSel, mEigval, mPvalues ] = Func_Sim_Bathia(n,d,NREP,Nboot)
% Essa função é utilizada para simular a estimação da dimensão de uma série
% temporal de curvas seguindo como fizeram Bathia et al (2010). As entradas
% são o tamanho amostral (número de curvas observadas) n, a dimensão do
% processo d, o número de réplicas da simulação NREP, o número de réplicas
% bootstrap Nboot no teste da dimensão.

% Parâmetros já definidos da simulação
p = 5;  % lag máximo utilizado
alpha = .05;  % nível de significância do teste bootstrap

% Outputs da função
mPvalues = ones(NREP,d);  % valores-p dos testes dos 8 maiores autovalores
mEigval = zeros(NREP,8);  % 8 maiores autovalores de D
vDimSel = (d+1)*ones(NREP,1);  % dimensões selecionadas em cada réplica


% processos para gerar a dependência temporal da função não observada
if d==2
    model_xi1 = arima('Constant',0,'AR',{-.65},'Variance',1.5);model_xi2 = arima('Constant',0,'AR',{.4},'Variance',1.5);
else
    if d==4
        model_xi1 = arima('Constant',0,'AR',{-.775},'Variance',1.5);model_xi2 = arima('Constant',0,'AR',{.65},'Variance',1.5);model_xi3 = arima('Constant',0,'AR',{-.525},'Variance',1.5);model_xi4 = arima('Constant',0,'AR',{.4},'Variance',1.5);
    else
        % deve ser o caso d==6
        model_xi1 = arima('Constant',0,'AR',{-.8167},'Variance',1.5);model_xi2 = arima('Constant',0,'AR',{.7333},'Variance',1.5);model_xi3 = arima('Constant',0,'AR',{-.65},'Variance',1.5);model_xi4 = arima('Constant',0,'AR',{.5667},'Variance',1.5);model_xi5 = arima('Constant',0,'AR',{-.4833},'Variance',1.5);model_xi6 = arima('Constant',0,'AR',{.4},'Variance',1.5);
    end
end

for kk=1:NREP
    
    % funcional observado
    fh_Y = cell(n,1);
    
    if d==2
        xi1 = simulate(model_xi1,n);xi2 = simulate(model_xi2,n);
        % gerando os coeficientes do funcional não observado
        for ii=1:n
            vZ = normrnd(0,1,10,1);
            fh_Y{ii} = @(x) xi1(ii)*sqrt(2)*cos(pi*x) + xi2(ii)*sqrt(2)*cos(2*pi*x)...
                + sqrt(2)*(vZ(1)*sin(pi*x) + vZ(2)*sin(2*pi*x)/2 ...
                + vZ(3)*sin(3*pi*x)/4 + vZ(4)*sin(4*pi*x)/8 ...
                + vZ(5)*sin(5*pi*x)/16 + vZ(6)*sin(6*pi*x)/32 ...
                + vZ(7)*sin(7*pi*x)/64 + vZ(8)*sin(8*pi*x)/128 ...
                + vZ(9)*sin(9*pi*x)/256 + vZ(10)*sin(10*pi*x))/512;
        end
    else
        if d==4
            xi1 = simulate(model_xi1,n);xi2 = simulate(model_xi2,n);
            xi3 = simulate(model_xi3,n);xi4 = simulate(model_xi4,n);
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
        else
            xi1 = simulate(model_xi1,n);xi2 = simulate(model_xi2,n);
            xi3 = simulate(model_xi3,n);xi4 = simulate(model_xi4,n);
            xi5 = simulate(model_xi5,n);xi6 = simulate(model_xi6,n);
            % deve ser o caso d==6
            for ii=1:n
                vZ = normrnd(0,1,10,1);
                fh_Y{ii} = @(x) xi1(ii)*sqrt(2)*cos(pi*x) + xi2(ii)*sqrt(2)*cos(2*pi*x)...
                    + xi3(ii)*sqrt(2)*cos(3*pi*x) + xi4(ii)*sqrt(2)*cos(4*pi*x)...
                    + xi5(ii)*sqrt(2)*cos(5*pi*x) + xi6(ii)*sqrt(2)*cos(6*pi*x)...
                    + sqrt(2)*(vZ(1)*sin(pi*x) + vZ(2)*sin(2*pi*x)/2 ...
                    + vZ(3)*sin(3*pi*x)/4 + vZ(4)*sin(4*pi*x)/8 ...
                    + vZ(5)*sin(5*pi*x)/16 + vZ(6)*sin(6*pi*x)/32 ...
                    + vZ(7)*sin(7*pi*x)/64 + vZ(8)*sin(8*pi*x)/128 ...
                    + vZ(9)*sin(9*pi*x)/256 + vZ(10)*sin(10*pi*x))/512;
            end
        end
    end
    
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
    
    % após obter D, calculamos os autovetores e autovalores da mesma.
    % as colunas de B terão os autovetores de D, que são os coeficientes de
    % ondaletas das autofunções associadas aos autovalores na diagonal de L
    [EigVecK,EigValK] = eig(Kstar);
    % 8 maiores autovalores de D
    mEigval(kk,:) = diag(EigValK(1:8,1:8));
    
    % calculando as 8 primeiras autofunções do operador analisado
    fh_OrtNormEigFunc = FuncGramSchm( fh_Ycent, EigVecK, d);
    mEta = integral(@(x) FuncMatProdYcentEigFunc( fh_Ycent, fh_OrtNormEigFunc, x ),0,1,'ArrayValued',true);
    
    checkDimSel = 0;
    d0 = 1;
    while d0<=d
        fh_OrtNormEigFunc = FuncGramSchm( fh_Ycent, EigVecK, d0);
        d_boot = FuncBootFuncDim(fh_Y, fh_OrtNormEigFunc, mEta(:,1:d0), d0, p, Nboot);
        mPvalues(kk,d0) = sum(d_boot>EigValK(d0+1,d0+1))/Nboot;
        if (mPvalues(kk,d0)>alpha)&&(checkDimSel==0)
            vDimSel(kk) = d0;
            checkDimSel = 1;
        end
        d0 = d0 + 1;
    end
        
end


end


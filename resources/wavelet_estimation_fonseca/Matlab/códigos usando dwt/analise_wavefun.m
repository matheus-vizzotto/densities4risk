% a funcao wavefun retorna os valores de psi e phi avaliadas nos pontos
% em XVAL
[PHI,PSI,XVAL] = wavefun('db2',20);

XVAL(471900)    % ponto x=0.45
PHI(471900)     % phi(0.45)
PSI(471900)     % psi(0.45)

% phi(0.45) avaliada com a funcao phi_n
phi_n(XVAL(471900),20,'db2')
% psi(0.45) avaliada com a funcao psi_n
psi_n(XVAL(471900),20,'db2')

% os valores de psi_n e PSI diferem para x=0.45
% PSI tem esse valor quando avaliada em x=0.45-1=-.55
XVAL(471900)-1
psi_n(XVAL(471900)-1,20,'db2')

% obtendo valores para construir o grafico da funcao de ondaletas
% usando a funcao psi_n
XV = linspace(-.99,1.99,1000);
%PHI_mf = zeros(length(XV),1);
PSI_mf = zeros(length(XV),1);
for ii = 1:length(XV)
    %PHI_mf(ii) = phi_n(XV(ii),20,'db2');
    PSI_mf(ii) = psi_n(XV(ii),20,'db2');
end

% comparando os graficos de psi obtidas da wavefun e de psi_n, vemos
% que a primeira esta deslocada em uma unidade para direita
subplot(2,1,1); plot(PSI); title('Psi obtido da funcao wavefun')
subplot(2,1,2); plot(XV,PSI_mf); title('Psi obtido da funcao psin')

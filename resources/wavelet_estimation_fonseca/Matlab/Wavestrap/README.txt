DimEstWav.m
% Neste código simulamos uma amostra de uma série temporal de funções como
% descrito na estudo numérico de Bathia, para analisar as funções
% DimEst_boot e DimEst_wavestrap
Esse é o código onde é feito a maioria dos testes das funções.


function [ vd0p1_boot ] = DimEst_boot( Y, NREP, B, p, N, wname)
% Função para realizar o procedimento bootstrap de Bathia et al (2010) para
% estimar a dimensão de um processo.
% o vetor vd0_boot conterá as dimensões bootstrap (d0+1) estimadas quando
% o processo tem dimensão d0 (número de colunas de B)
% Y é a matriz de funcionais observados, nas colunas, para diferentes
% dias, nas linhas.
% NREP é o número de réplicas bootstrap. p é lag máximo. wname é o nome
% da base de ondaletas utilizada.


function [ vd0p1_boot ] = DimEst_boot( Y, NREP, B, p, N, wname)
% Função para realizar o procedimento bootstrap de Bathia et al (2010) para
% estimar a dimensão de um processo.
% o vetor vd0_boot conterá as dimensões bootstrap (d0+1) estimadas quando
% o processo tem dimensão d0 (número de colunas de B)
% Y é a matriz de funcionais observados, nas colunas, para diferentes
% dias, nas linhas.
% NREP é o número de réplicas bootstrap. p é lag máximo. wname é o nome
% da base de ondaletas utilizada.


function [ vd0p1_boot ] = DimEst_wavestrap( A, NREP, B, p)
% Nessa função é realizado o procedimento bootstrap para estimar a 
% dimensão do processo, mas diferentemente do que é feito em DimEst_boot.m,
% aqui aproveitamos as decomposições obtidas para realizar a reamostragem
% nos termos formados com os coeficientes, o que deixa o código mais rápido.
% A é a matriz de coeficientes de ondaletas dos funcionais observados, 
% nas colunas, para diferentes dias, nas linhas.
% NREP é o número de réplicas bootstrap. p é lag máximo. 


function [ vd0p1_boot ] = DimEst_wavestrap2( Y, Nboot, B, p, N, wname)
% Nessa função é realizado um procedimento bootstrap para estimar a 
% dimensão do processo. Obtemos os coeficientes de ondaleta para Y e
% realizamos uma limiarização rígida, identificando quais destes coef. são
% ruído. A reamostragem é realizada sobre os coef. limiarizados
% (padronizados anteriormete), que são sorteados para substituir os zeros
% da matriz de coef. limiarizados, nos respectivos dias e níveis.
% Y é a matriz de funcionais observados, B é a matriz cujas linhas são
% observações das funções que geram Y. N e wname são argumentos da função
% waavdec. Nboot é o número de réplicas bootstrap. p é lag máximo. 


function [ vd0p1_boot ] = DimEst_wavestrap3( A, Lw, Nboot, B, p, N, wname)
% Nessa função é realizado o procedimento bootstrap para estimar a 
% dimensão do processo, as decomposições obtidas são aproveitadas na
% reamostragem, que teoricamente é feita sobre o resíduo obtido após ser 
% feita uma limiarização sobre os coeficientes das decomposições da média
% do processo e das autofunções.
% A é a matriz de coeficientes de ondaletas dos funcionais observados, 
% nas colunas, para diferentes dias, nas linhas.
% NREP é o número de réplicas bootstrap. p é lag máximo. 


function [ vd0p1_boot ] = DimEst_wavestrap_thresh( A, NREP, B, p, N, wname)
% Nessa função é realizado o procedimento bootstrap para estimar a 
% dimensão do processo, as decomposições obtidas são aproveitadas na
% reamostragem, que teoricamente é feita sobre o resíduo obtido após ser 
% feita uma limiarização sobre os coeficientes das decomposições da média
% do processo e das autofunções.
% A é a matriz de coeficientes de ondaletas dos funcionais observados, 
% nas colunas, para diferentes dias, nas linhas.
% NREP é o número de réplicas bootstrap. p é lag máximo. 


function [Aboot] = ResampWavestrap(A, Lw, Nboot, N, wname)
% função para gerar os coeficientes de ondaletas de uma amostra
% bootstrap que foi gerada através do processo que usa limiarização
% para identificar quais coeficientes são ruído ou não na amostra
% original. O argumento A é a matriz dos coeficientes da amostra
% original, Lw é retornado da função wavedec enquanto N e wname
% são argumentos de wavedec.


gerar_amostra_via_WaveCoeff.m
% Nesse código é feito o esquema para gerar amostras bootstrap do funcional
% observado Y, usando limiarização para identificar quais coeficientes são
% ruídos, e realizando a reamostragem nestes.
Esse esquema de reamostragem está implementado na função DimEst_wavestrap2


comparar_IC_d_WaveBoot.m
% Neste código simulamos uma amostra de uma série temporal de funções como
% descrito na estudo numérico de Bathia. Depois, geramos versões bootstrap
% da amostra observada Y usando o seguinte procedimento: os coeficientes
% das decomposições de cada dia são limiarizados com hard thresholding, onde
% identificamos os coef. que são ruído; esses coeficientes são repadronizados
% em cada nível para terem média zero e desvio padrão um; a reamostragem é
% feita sobre esses rezíduos padronizados, que são reescalados pela média e std
% dos níveis e colocados no lugar do coef. limiarizados. Com essas amostras
% bootstrap obtemos vetores de autovalores calculados, e então contruímos
% intervalos de confianças, que tem seus comprimentos comparados


gerar_amostra_via_WaveCoeff.m
% Nesse código é feito o esquema para gerar amostras bootstrap do funcional
% observado Y, usando limiarização para identificar quais coeficientes são
% ruídos, e realizando a reamostragem nestes.




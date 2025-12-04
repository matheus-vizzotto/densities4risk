% Arquivo para ler e analisar os dados simulados

fn100d2 = matfile('file_Sim_DimEstWav_n100d2.mat');
fn300d2 = matfile('file_Sim_DimEstWav_n300d2.mat');
fn600d2 = matfile('file_Sim_DimEstWav_n600d2.mat');

fn100d4 = matfile('file_Sim_DimEstWav_n100d4.mat');
fn300d4 = matfile('file_Sim_DimEstWav_n300d4.mat');
fn600d4 = matfile('file_Sim_DimEstWav_n600d4.mat');

fn100d6 = matfile('file_Sim_DimEstWav_n100d6.mat');
fn300d6 = matfile('file_Sim_DimEstWav_n300d6.mat');
fn600d6 = matfile('file_Sim_DimEstWav_n600d6.mat');

% Gráfico com a porcentagem das dimensões selecionadas para cada tamanho 
% amostral e  valor verdadeiro de d

mTbn100d2 = tabulate(fn100d2.vDimSel);mTbn300d2 = tabulate(fn300d2.vDimSel);mTbn600d2 = tabulate(fn600d2.vDimSel);
mTbn100d4 = tabulate(fn100d4.vDimSel);mTbn300d4 = tabulate(fn300d4.vDimSel);mTbn600d4 = tabulate(fn600d4.vDimSel);
mTbn100d6 = tabulate(fn100d6.vDimSel);mTbn300d6 = tabulate(fn300d6.vDimSel);mTbn600d6 = tabulate(fn600d6.vDimSel);

subplot(3,1,1)
bar([mTbn100d2(1:3,3)';mTbn300d2(1:3,3)';mTbn600d2(1:3,3)'])
legend('1','2','3','Location','northwest');legend('boxoff');
title('(a) d=2');
xt = {'n=100';'n=300';'n=600'};
set(gca,'xticklabel',xt, 'TickLabelInterpreter', 'latex');
subplot(3,1,2)
bar([mTbn100d4(3:5,3)';mTbn300d4(3:5,3)';mTbn600d4(3:5,3)'])
legend('3','4','5','Location','northwest');legend('boxoff');
title('(b) d=4');
xt = {'n=100';'n=300';'n=600'};
set(gca,'xticklabel',xt, 'TickLabelInterpreter', 'latex');
subplot(3,1,3)
bar([mTbn100d6(5:7,3)';mTbn300d6(5:7,3)';mTbn600d6(5:7,3)'])
legend('5','6','7','Location','northwest');legend('boxoff');
title('(c) d=6');
xt = {'n=100';'n=300';'n=600'};
set(gca,'xticklabel',xt, 'TickLabelInterpreter', 'latex');

% Gráfico dos 8 maiores autovalores para cada tamanho amostral e valor
% verdadeiro de d

subplot(3,1,1)
plot(1:8,mean(fn100d2.mEigval)/(10^5),'-.r*',...
1:8,mean(fn300d2.mEigval)/(10^5),'--mo',1:8,mean(fn600d2.mEigval)/(10^5),':bs')
title('(a) d=2');
ylim([0 10]);xlim([1 8]);
legend('n=100','n=300','n=600');legend('boxoff');
subplot(3,1,2)
plot(1:8,mean(fn100d4.mEigval)/(10^5),'-.r*',...
1:8,mean(fn300d4.mEigval)/(10^5),'--mo',1:8,mean(fn600d4.mEigval)/(10^5),':bs')
title('(b) d=4');
ylim([0 10]);xlim([1 8]);
legend('n=100','n=300','n=600');legend('boxoff');
subplot(3,1,3)
plot(1:8,mean(fn100d6.mEigval)/(10^5),'-.r*',...
1:8,mean(fn300d6.mEigval)/(10^5),'--mo',1:8,mean(fn600d6.mEigval)/(10^5),':bs')
title('(c) d=6');
ylim([0 10]);xlim([1 8]);
legend('n=100','n=300','n=600');legend('boxoff');

% boxplot dos valores-p dos testes do d-ésimo e (d+1)-ésimo autovalor, para
% cada tamanho amostral e valor verdadeiro de d

mPvn100d2 = fn100d2.mPvalues;mPvn300d2 = fn300d2.mPvalues;mPvn600d2 = fn600d2.mPvalues;
mPvn100d4 = fn100d4.mPvalues;mPvn300d4 = fn300d4.mPvalues;mPvn600d4 = fn600d4.mPvalues;
mPvn100d6 = fn100d6.mPvalues;mPvn300d6 = fn300d6.mPvalues;mPvn600d6 = fn600d6.mPvalues;

subplot(2,1,1)
aux = [mPvn100d2(:,1),mPvn300d2(:,1),mPvn600d2(:,1),...
    mPvn100d4(:,3),mPvn300d4(:,3),mPvn600d4(:,3),...
    mPvn100d6(:,5),mPvn300d6(:,5),mPvn600d6(:,5)];
xt = {'\begin{tabular}{c}d=2\\n=100\end{tabular}';'\begin{tabular}{c}d=2\\n=300\end{tabular}';'\begin{tabular}{c}d=2\\n=600\end{tabular}';...
    '\begin{tabular}{c}d=4\\n=100\end{tabular}';'\begin{tabular}{c}d=4\\n=300\end{tabular}';'\begin{tabular}{c}d=4\\n=600\end{tabular}';...
    '\begin{tabular}{c}d=6\\n=100\end{tabular}';'\begin{tabular}{c}d=6\\n=300\end{tabular}';'\begin{tabular}{c}d=6\\n=600\end{tabular}'};
boxplot(aux);
line([0:10],.1*ones(11,1),'LineStyle','--');
set(gca,'xticklabel',xt, 'TickLabelInterpreter', 'latex');
title('(a) Valores-p para o d-ésimo maior autovalor');
subplot(2,1,2)
aux = [mPvn100d2(:,2),mPvn300d2(:,2),mPvn600d2(:,2),...
    mPvn100d4(:,4),mPvn300d4(:,4),mPvn600d4(:,4),...
    mPvn100d6(:,6),mPvn300d6(:,6),mPvn600d6(:,6)];
xt = {'\begin{tabular}{c}d=2\\n=100\end{tabular}';'\begin{tabular}{c}d=2\\n=300\end{tabular}';'\begin{tabular}{c}d=2\\n=600\end{tabular}';...
    '\begin{tabular}{c}d=4\\n=100\end{tabular}';'\begin{tabular}{c}d=4\\n=300\end{tabular}';'\begin{tabular}{c}d=4\\n=600\end{tabular}';...
    '\begin{tabular}{c}d=6\\n=100\end{tabular}';'\begin{tabular}{c}d=6\\n=300\end{tabular}';'\begin{tabular}{c}d=6\\n=600\end{tabular}'};
boxplot(aux);
line([0:10],.1*ones(11,1),'LineStyle','--');
set(gca,'xticklabel',xt, 'TickLabelInterpreter', 'latex');
title('(b) Valores-p para o (d+1)-ésimo maior autovalor');


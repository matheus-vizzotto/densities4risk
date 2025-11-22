function psiv = psi_n(x,n,wname)
% função para calcular a função de ondaleta psi avaliada em x com aproximação
% de ordem n nos diáticos e para um filtro de ondaletas 'wname'
h = fliplr(wfilters(wname));
N = length(h)/2;

if(x<=1-N || x>=N)
    psiv = 0;
else
    psiv = 0;
    for ii = (-1):(2*N-2)
        psiv = psiv + ((-1)^ii)*sqrt(2)*h(2+ii)*phi_n(2*x+ii,n,wname);
    end
end

end


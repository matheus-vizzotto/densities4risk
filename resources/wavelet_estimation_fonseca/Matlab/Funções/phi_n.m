function phiv = phi_n(x,n,wname)
% função para calcular a função escala avaliada em x com aproximação
% de ordem n nos diáticos e para um filtro de ondaletas 'wname'

h = fliplr(wfilters(wname));
N = length(h)/2;
T0 = zeros(2*N-1,2*N-1);
T1 = zeros(2*N-1,2*N-1);

for ii = 1:(2*N-1)
for jj = 1:(2*N-1)
	if (2*ii-jj-1<0) h_ind=0; elseif (2*ii-jj-1>2*N-1) h_ind=0;... 
    else h_ind=h(2*ii-jj); end
	T0(ii,jj) = sqrt(2)*h_ind;
	if (2*ii-jj<0) h_ind=0;elseif (2*ii-jj>2*N-1) h_ind=0;...
    else h_ind=h(2*ii-jj+1); end
	T1(ii,jj) = sqrt(2)*h_ind;
end
end

  if (x<=0 || x>=length(h)-1)
    phiv = 0;
  else
    k = 1;
    y = x;
    while (y>=0)
      if (y<1)
          vd = dyadn(y,n);
          if (vd(1)==0) mT = T0; else mT = T1; end
          for ii = 2:length(vd)
              if (vd(ii)==0) mT = mT*T0; else mT = mT*T1; end
          end
          phiv = mT(k,1);
          break
      else
          y = y - 1;
          k = k + 1;
      end
    end
  end
end

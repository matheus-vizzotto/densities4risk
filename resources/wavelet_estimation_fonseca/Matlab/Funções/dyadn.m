function d = dyadn(x,n)
% retorna os n primeiros números da expansão diádica
% de um número x em [0,1]
  d = ones(n,1);
  xd = x;
  if xd==1
      d = d*9;
  else
  for ii = 1:n
    d(ii) = floor(2*xd);
    xd = 2*xd - floor(2*xd);
  end
  end
end

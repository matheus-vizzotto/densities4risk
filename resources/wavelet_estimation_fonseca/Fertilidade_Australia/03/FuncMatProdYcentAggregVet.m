function [ y ] = FuncMatProdYcentAggregVet( C, w, p, k)


% tamanho da amostra
[nt,n] = size(C);
% número de pesos no vetor w
delta = length(w);

% obtendo os funcionais observados ponderados
mY = zeros(delta,delta,nt);
mZ = zeros(n-p-delta+1,nt);
for jj = 1:(n-p-delta+1)
    for ii = 0:(delta-1)
        mY(ii+1,:,:) = C(:,(jj+k-ii):(jj+k+delta-1-ii))';
        %for ll=0:(delta-1)
        %    mY(ii+1,ll+1,:) = C(:,jj+k+ll-ii);
        %end
    end
    if delta==1
        mZ(jj,:) = (w^2)*mY(1,1,:);
    else
            for ii = 1:nt
                mZ(jj,ii) = w'*mY(:,:,ii)*w;
            end
    end
end

% matriz com produtos dos funcionais observados ponderados
y = mZ*mZ';

end


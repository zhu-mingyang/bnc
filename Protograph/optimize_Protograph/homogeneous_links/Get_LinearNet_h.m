function [h] = Get_LinearNet_h(k,e,M,q,zeta_table)
% k: the length of linear network
% e: erasure rate vector of length k
% M: batch size
% q: the number of elements in the base field
h_tmp = zeros(k,M+1);
for r = 0:M
    h_tmp(1,r+1) = nchoosek(M,r)*(1-e(1))^r*e(1)^(M-r);
end
for it = 2:k
    for r = 0:M
        for i = r:M
            for j = r:M
                h_tmp(it,r+1) = h_tmp(it,r+1) + h_tmp(it-1,i+1)*nchoosek(M,j)*(1-e(it))^j*e(it)^(M-j)*zeta(i,j,r,q,zeta_table);
            end
        end
    end
end
h = h_tmp(k,:);
end

function [y] = zeta(d,k,r,q,zeta_table)  % zeta_r^{d,k}
y = zeta_table(d+1,r+1)*zeta_table(k+1,r+1)/zeta_table(r+1,r+1)/q^((d-r)*(k-r));
end
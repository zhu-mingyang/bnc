function [P,full_rk] = GaussianEliminate(H,q)
% H must be a tall matrix, i.e., size(H,1) >= size(H,2)
% if H has full rank, P is the transfer matrix, full_rk = 1
% if H does not have full rank, P is invalid, full_rk = 0
[M,N] = size(H);
P = gf(eye(M),log2(q));
for i = 1:N
    rp = M - i + 1;
    for j = 1:M - i + 1
        if (H(rp, N-i+1) ~= 0)
            break;
        else
            rp = rp - 1;
        end
    end
    if (rp == 0)                                                            % the tall matrix is not of full rank
       P = [];
       full_rk = 0;
       return;
    end
    
    if (rp ~= M - i + 1)
        H(M - i + 1,:) =  H(M - i + 1,:)+H(rp,:);
        P(M - i + 1,:) =  P(M - i + 1,:)+P(rp,:);
    end
    
    for j = 1:M
       if (j ~= M - i + 1 && H(j,N-i+1) ~= 0) 
           if (H(M-i+1,N-i+1) == 0)
               error('ERROR: Gaussian Elimination!')
           end
           factor = H(M-i+1,N-i+1)^(-1)*H(j,N-i+1);
           H(j,:) = H(j,:)+factor*H(M-i+1,:);
           P(j,:) = P(j,:)+factor*P(M-i+1,:);
       end
    end
    factor = H(M-i+1,N-i+1)^(-1);
    H(M-i+1,:) = factor*H(M-i+1,:);
    P(M-i+1,:) = factor*P(M-i+1,:);
end
full_rk = 1;
end
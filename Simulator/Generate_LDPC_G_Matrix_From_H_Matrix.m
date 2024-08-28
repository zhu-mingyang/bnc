function [H_reg,H_sys,G_sys,pos] = Generate_LDPC_G_Matrix_From_H_Matrix(H_gf,q,m,n)
% *****************REVISION Begin 2024/06/24******************************
% pos needs to be output. If the columns of the precode H matrix have been adjusted, 
% we need to adjust the batch profile (otherwise the performance will slightly degrade)
% *****************REVISION END******************************

k = n-m;
[H_sys,H_reg,~,pos] = GE(H_gf);
if (m > 0)
    G_sys = [gf(eye(k),log2(q)) H_sys(:,1:k)'];
else
    G_sys = gf(eye(k),log2(q));
end
end

function [H_sys, H_new, fs, pos] = GE(H_input)
fs = 0;
[M,N] = size(H_input);
H_sys = H_input;
pos = 1:N;
for i = 1:M
    rp = M - i + 1;
    for j = 1:M - i + 1
        if (H_sys(rp, N-i+1) ~= 0)
            break;
        else
            rp = rp - 1;
        end
    end
    if (rp == 0)
       [rp, scp] = SeekColumn(H_sys, i, M, N);
       if (rp == 0)
           error('H matrix is not of full rank!');
       end
       H_sys = SwitchColumn(H_sys, N-i+1, scp);
       fs = fs + 1;
       temp = pos(N-i+1);
       pos(N-i+1) = pos(scp);
       pos(scp) = temp;
    end
    
    if (rp ~= M - i + 1)
        H_sys(M - i + 1,:) =  H_sys(M - i + 1,:) + H_sys(rp,:);
    end
    
    for j = 1:M
       if (j ~= M - i + 1 && H_sys(j,N-i+1) ~= 0) 
           if (H_sys(M-i+1,N-i+1) == 0)
               error('ERROR: Gaussian Elimination!')
           end
           factor = H_sys(M-i+1,N-i+1)^(-1)*H_sys(j,N-i+1);
           H_sys(j,:) = H_sys(j,:)+factor*H_sys(M-i+1,:);
       end
    end
    H_sys(M-i+1,:) = H_sys(M-i+1,N-i+1)^(-1)*H_sys(M-i+1,:);
end
H_new = H_input(:,pos);
end

function [rp, scp] = SeekColumn(H, cp, M, N)
rp = 0;
for i = N - cp:-1:1
   for j = M - cp + 1:-1:1
      if (H(j,i) ~= 0)
          rp = j;
          scp = i;
          return;
      end
   end
end
end

function [H] = SwitchColumn(H, cp, scp)
    temp = H(:,cp);
    H(:,cp) = H(:,scp);
    H(:,scp) = temp;
end



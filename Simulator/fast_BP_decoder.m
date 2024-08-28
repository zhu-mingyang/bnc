function [num_undecoded] = fast_BP_decoder(A,G,H,q,K,n,m)
% joint BP decoder for precoded BATS codes with n standard batches and m precoded batches

decoded_VN = zeros(1,K);                                                  % the decoder exits when sum(decoded_flag)==K

GH = cell(1,n);     
rk = zeros(1,m+n);     % main revision: store the rank information, avoid computing rank() many times in decoding
for i = 1:n
    GH{i} = G{i}*H{i};
    rk(i) = rank(GH{i});
end
for i = n+1:n+m
    rk(i) = 1;    
end

while (1)
    if (sum(decoded_VN)==K)
        break;
    end
    decodable = 0;                                                          % indicate if we can find a decodable batch
    G_update_pos = zeros(1,n+m);
    for i = 1:n+m
        if (isempty(A{i})==1)                                              % this batch has been decoded
            continue;
        end   
        dg = max(size(A{i}));                
        if (rk(i) == dg)
            if (i <= n)     % for standard batches
                A_ = GH{i}';
                if (dg~=size(A_,2))
                    error('matrix size error.')
                end
                %b = Y{i}';
                %[P,~] = GaussianEliminate(A_,q);                                % solve (GH)'B'=Y'
                decoded_VN(A{i}) = 1;
                %B_temp = P*b;
                %B_temp = B_temp(M-dg+1:M,:)';                                   % our GaussianEliminate function puts I at the bottom
                %B_dec(:,A{i}) = B_temp;
            else            % for precoded batches
                decoded_VN(A{i}) = 1;
                %B_temp = Y{i}*(A_^(-1));
                %B_dec(:,A{i}) = Y{i}*(A_^(-1));
            end
            for j = 1:n+m
                if (isempty(A{j})==1 || j==i)
                    continue;
                end
                for k = 1:dg
                    s = find(A{j}==A{i}(k));
                    if (isempty(s)==0)
                        A{j}(:,s) = [];                        
                        %Y{j} = Y{j}-B_temp(:,k)*G{j}(s,:)*H{j};
                        G_temp = G{j}.x;
                        G_temp(s,:) = [];                                   % it seems that MATLAB cannot delete row/column of gf matrix
                        G{j} = gf(G_temp,log2(q));
                        G_update_pos(j) = 1;
                    end
                end
            end
            A{i} = [];
            rk(i) = 0;  % has been decoded
            decodable = 1;
            for j = 1:n  % update rank information for standard batches
                if (G_update_pos(j)==1 && isempty(A{j})==0)                   
                    GH{j} = G{j}*H{j};
                    rk(j) = rank(GH{j});
                end
            end
            break;
        end
    end
    if (decodable==0)
        break;
    end
end
num_undecoded = K-sum(decoded_VN);
end


% function [num_undecoded] = fast_BP_decoder(A,G,H,q,K,n,m)   % copy from V5, slower than the new version
% % joint BP decoder for precoded BATS codes with n standard batches and m precoded batches
% decoded_VN = zeros(1,K);                                                  % the decoder exits when sum(decoded_flag)==K
% 
% while (1)
%     if (sum(decoded_VN)==K)
%         break;
%     end
%     decodable = 0;                                                          % indicate if we can find a decodable batch
%     for i = 1:n+m
%         if (isempty(A{i})==1)                                              % this batch has been decoded
%             continue;
%         end
%         A_ = G{i}*H{i};        
%         dg = max(size(A{i}));                
%         if (rank(A_) == dg)
%             if (i <= n)     % for standard batches
%                 A_ = A_';
%                 if (dg~=size(A_,2))
%                     error('matrix size error.')
%                 end
%                 %b = Y{i}';
%                 %[P,~] = GaussianEliminate(A_,q);                                % solve (GH)'B'=Y'
%                 decoded_VN(A{i}) = 1;
%                 %B_temp = P*b;
%                 %B_temp = B_temp(M-dg+1:M,:)';                                   % our GaussianEliminate function puts I at the bottom
%                 %B_dec(:,A{i}) = B_temp;
%             else            % for precoded batches
%                 decoded_VN(A{i}) = 1;
%                 %B_temp = Y{i}*(A_^(-1));
%                 %B_dec(:,A{i}) = Y{i}*(A_^(-1));
%             end
%             for j = 1:n+m
%                 if (isempty(A{j})==1 || j==i)
%                     continue;
%                 end
%                 for k = 1:dg
%                     s = find(A{j}==A{i}(k));
%                     if (isempty(s)==0)
%                         A{j}(:,s) = [];                        
%                         %Y{j} = Y{j}-B_temp(:,k)*G{j}(s,:)*H{j};
%                         G_temp = G{j}.x;
%                         G_temp(s,:) = [];                                   % it seems that MATLAB cannot delete row/column of gf matrix
%                         G{j} = gf(G_temp,log2(q));
%                     end
%                 end
%             end
%             A{i} = [];
%             decodable = 1;
%             break;
%         end
%     end
%     if (decodable==0)
%         break;
%     end
% end
% num_undecoded = K-sum(decoded_VN);
% end
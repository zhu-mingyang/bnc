function [B_dec,fail,INAC_num] = BP_Decoder(A,G,H,Y,q,K,n,m,M,N_INAC)
% Packet size is fixed to 1, i.e., one symbol
% joint BP decoder for precoded BATS codes with n standard batches and m precoded batches
decoded_VN = zeros(1,K);                                                  % the decoder exits when sum(decoded_flag)==K
B_dec = gf(zeros(1,K),log2(q));
B_dec_with_ina = gf(zeros(N_INAC+1,K),log2(q));     % represent decoded symbols by inactive symbols, every column denotes a decoded symbol = a0+a1b1+a2b2+...
INAC_sym_set = zeros(1,K);          % the index of inactive symbols
INAC_num = 0;                       % number of inactive symbols actually used
INAC_lcons_num = 0;            % number of linear constraints for inactive symbols
INAC_lcons = gf(zeros(N_INAC+1,2*N_INAC),log2(q)); % linear constraint matrix, every column is a linear constraint on a0+a1b1+a2b2+...=0
while (1)
    if (sum(decoded_VN)==K)
        break;
    end
    decodable = 0;                                                          % indicate if we can find a decodable batch
    for i = 1:n+m
        if (isempty(A{i})==1)                                              % this batch has been decoded
            continue;
        end
        A_ = G{i}*H{i};        
        dg = max(size(A{i}));                
        if (rank(A_) == dg)
            if (i <= n)     % for standard batches
                A_ = A_';
                if (dg~=size(A_,2))
                    error('matrix size error.')
                end
                b = Y{i}';
                [P,~] = GaussianEliminate(A_,q);                                % solve (GH)'B'=Y'
                decoded_VN(A{i}) = 1;
                B_temp = P*b;
                B_temp = B_temp(M-dg+1:M,:)';                                   % our GaussianEliminate function puts I at the bottom
                B_dec_with_ina(:,A{i}) = B_temp;
            else            % for precoded batches
                decoded_VN(A{i}) = 1;
                B_temp = Y{i}*(A_^(-1));
                B_dec_with_ina(:,A{i}) = Y{i}*(A_^(-1));
            end
            for j = 1:n+m
                if (isempty(A{j})==1 || j==i)
                    continue;
                end
                for k = 1:dg
                    s = find(A{j}==A{i}(k));
                    if (isempty(s)==0)
                        A{j}(:,s) = [];                        
                        Y{j} = Y{j}-B_temp(:,k)*G{j}(s,:)*H{j};
                        G_temp = G{j}.x;
                        G_temp(s,:) = [];                                   % it seems that MATLAB cannot delete row/column of gf matrix
                        G{j} = gf(G_temp,log2(q));
                    end
                end
            end
            A{i} = [];
            Y{i} = Y{i} - B_temp*G{i}*H{i};
            decodable = 1;
            break;
        end
    end
    if (decodable==0)
        if (INAC_num < N_INAC)          % start inactivation decoding
            for i = 1:K
                if (decoded_VN(i) == 0)
                    decoded_VN(i) = 1;          % inactivate this symbol and set to decoded
                    INAC_num = INAC_num + 1;
                    INAC_sym_set(INAC_num) = i;
                    B_dec_with_ina(INAC_num+1,i) = 1;
                    for j = 1:m+n
                        if (isempty(A{j})==1)   % the batch was deleted
                            continue;
                        end
                        s = find(A{j}==i);      % find a batch containing the inactive symbol
                        if (isempty(s)==0)
                            A{j}(:,s) = [];                        
                            Y{j}(INAC_num+1,:) = G{j}(s,:)*H{j};
                            G_temp = G{j}.x;
                            G_temp(s,:) = [];                                   % it seems that MATLAB cannot delete row/column of gf matrix
                            G{j} = gf(G_temp,log2(q));
                        end
                    end
                    break;
                end
            end
        else
            break;
        end
    end
end
if (sum(decoded_VN)==K)    
    if (INAC_num == 0)
        fail = 0;
        B_dec = B_dec_with_ina(1,:);
    else            % solve the linear constaints for inactive symbols
        for i = 1:m+n   % find constraints for inactive symbols
            for j = 1:size(Y{i},2)
                if (sum(Y{i}(:,j)~=0)>0)
                    INAC_lcons_num = INAC_lcons_num + 1;
                    INAC_lcons(:,INAC_lcons_num) = Y{i}(:,j);
                end
            end
        end
        A_ = INAC_lcons(2:INAC_num+1,1:INAC_lcons_num)';
        if (rank(A_) == INAC_num)
            fail = 0;
            b = INAC_lcons(1,1:INAC_lcons_num)';
            [P,~] = GaussianEliminate(A_,q);
            tmp = P*b;
            base = [gf(1,log2(q));tmp(INAC_lcons_num-INAC_num+1:INAC_lcons_num)];     % (1,b1,b2,...)
            for i = 1:K
                for j = 1:INAC_num+1
                    B_dec(i) = B_dec(i) + B_dec_with_ina(j,i)*base(j);
                end
            end
        else    % inactivation decoding fails
            fail = 1;
            B_dec = B_dec_with_ina(1,:);
        end
    end
else    % BP decoding fails
    fail = 1;       
    B_dec = B_dec_with_ina(1,:);
end
end
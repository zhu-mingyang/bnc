%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the program for optimization of protograph based on Homogeneous Links
% The code is based on the paper: M. Zhu et al., Protograph-Based Batched Network Codes

% NOTICE: The optimization method used in the program is very simple, primarily based on random and brute-force searches. 
%         It is merely a straightforward application of the theory proposed in the paper. By referring to the optimization methods 
%         for LDPC protographs, it is not difficult to find more efficient optimization approaches.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc;

%%%%%%%%%%%%%%%%%%%%%%%%% Variables that need to be adjusted %%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = 8;         % batch size
q = 256;        % base field
n_hops = 3;     % the length of line network (i.e., the number of hops)

% protograph of the precode (B^(1), this is fixed in the optimization)
B_p = [ 1 3 1 1 1 1 1 0;
        1 3 2 0 1 0 0 1;
        0 1 2 1 1 1 1 1];
% B_p = [1 3 2 1 1 1 2 1; 
%        1 3 2 1 1 1 1 2];

n_v = 8;        % number of variable nodes
n1 = 6;         % number of rows in the core protomatrix (B_c^(2))
n2 = 6;         % number of rows in the extension protomatrix (B_e^(2))

delta_c = [0.74 0.86 0.86 0.86 0.92 0.92];      % the puncturing vector for the core protomatrix (this will be optimized)
delta_e = [0.88 0.88 0.8 0.8 0.8 0.8];          % the puncturing vector for the extension protomatrix (this is fixed)

row_core_deg = [6 8 10 13 15 19];             % initial row degrees of the core protomatrix


% Following is an empirical approach which has not been mentioned in the paper
% This approach assigns different probability distributions to different elements of the matrix. 
% This may speed up the traversal process.
light_candidate = [0,1,2];
light_prob = [0.25 0.45 0.3];
dark_candidate = [0,1,2,3,4];
dark_prob = [0 0.2 0.2 0.3 0.3];
dark_area = [0 1 1 0 0 0 1 1;
             0 1 1 0 0 0 1 1;
             0 1 1 0 0 0 1 1;
             0 1 1 0 0 0 1 1;
             1 1 1 1 1 1 1 1;
             1 1 1 1 1 1 1 1];
extension_candidate = [0,1,2];
extension_prob = [1/3 1/3 1/3];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End of Variables that need to be adjusted %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialize some variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zeta = generate_zeta(M,q);          % store the values of zeta function
th_opt = 0;                      % In this program, the Threshold is measured by erasure probability, so initialize it to 0
B_bats_opt = gen_random_protograph(n1,n_v,row_core_deg,light_candidate,dark_candidate,light_prob,dark_prob,dark_area);
punc_vector_opt = delta_c;
punc_sum = sum(delta_c);


%%%%%%%%%%%%%%%%%%%%%%%%%%%% Iterative Optimization of Core Matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iter = 1:100
    find_a_new_protograph = 0;
    row_core_deg = sum(B_bats_opt,2);
    for iter_proto = 1:10       % fix row degrees to optimize
        th = zeros(1,50);
        parfor i = 1:50
            B_bats = gen_random_protograph(n1,n_v,row_core_deg,light_candidate,dark_candidate,light_prob,dark_prob,dark_area);
            B_record(i,:) = reshape(B_bats,1,n1*n_v);
            th(i) = Proto_BATS_Threshold_Approximate(B_bats,B_p,M,q,zeta,punc_vector_opt,n_hops);
        end
        [th_max,th_max_id] = max(th);
        if (th_max > th_opt)
            find_a_new_protograph = 1;
            th_opt = th_max;
            B_bats_opt = reshape(B_record(th_max_id,:),n1,n_v);            
        end
    end
    col_core_deg = sum(B_bats_opt);
    for iter_proto = 1:10       % fix column degrees to optimize
        th = zeros(1,50);
        parfor i = 1:50
            B_bats = gen_random_protograph(n_v,n1,col_core_deg,light_candidate,dark_candidate,light_prob,dark_prob,dark_area');
            B_bats = B_bats';           
            B_record(i,:) = reshape(B_bats,1,n1*n_v);

            % The following condition is empirical: For finite-length BATS codes, if no degree < M, the performance may be poor
            if (sum(sum(B_bats,2) < M) > 0)
                th(i) = Proto_BATS_Threshold_Approximate(B_bats,B_p,M,q,zeta,punc_vector_opt,n_hops);
            else
                th(i) = -1;
            end
        end
        [th_max,th_max_id] = max(th);
        if (th_max > th_opt)
            find_a_new_protograph = 1;
            th_opt = th_max;
            B_bats_opt = reshape(B_record(th_max_id,:),n1,n_v);            
        end
    end
    if (find_a_new_protograph == 0)
        continue;
    end
    for iter_punc = 1:4
        punc_vector_opt_tmp = punc_vector_opt;
        for j = 1:n1        % adjust the elements of puncturing vector one by one  
            punc_vector = zeros(51,n1);
            num_valid_punc_vec = 0;
            th = zeros(1,51);
            delta_max = punc_vector_opt(j)*0.1;
            delta_step = 2*delta_max/50;
            for s = 1:51
                punc_vector_tmp = punc_vector_opt;
                delta = (s-1)*delta_step - delta_max;% [-delta_max, +delta_max]              
                punc_vector_tmp(j) = punc_vector_tmp(j) + delta;
                punc_vector_tmp = punc_vector_tmp/(delta+punc_sum)*punc_sum;
                if (sum(punc_vector_tmp > 1) == 0 && sum(punc_vector_tmp < 0) == 0)
                    num_valid_punc_vec = num_valid_punc_vec + 1;
                    punc_vector(num_valid_punc_vec,:) = punc_vector_tmp;
                end
            end
            if (num_valid_punc_vec == 0)
                error('num_valid_punc_vec = 0!');
            end
            parfor ss = 1:num_valid_punc_vec
                th(ss) = Proto_BATS_Threshold_Approximate(B_bats_opt,B_p,M,q,zeta,punc_vector(ss,:),n_hops);
            end
            [th_max,th_max_id] = max(th);
            if (th_max > th_opt)
                th_opt = th_max;
                punc_vector_opt_tmp = punc_vector(th_max_id,:);
            end
        end
        punc_vector_opt = punc_vector_opt_tmp;
    end
    rate = (size(B_bats_opt,2)-size(B_p,1))/(size(B_bats_opt,1)-sum(punc_vector_opt));
    th_capacity = capacity(n_hops,th_opt*ones(1,n_hops),M,q,zeta);
    fprintf('Core Protomatrix Optimization, iter = %d, th_erasure = %.4f, th_capacity = %.4f, Rate = %.4f\n',iter,th_opt,th_capacity,rate);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Iterative Optimization of Extension Matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i_ext = 1:n2
    punc_vector_opt = [punc_vector_opt,delta_e(i_ext)];
    B_bats_opt_old = B_bats_opt;
    for iter = 1:100
        find_a_new_protograph = 0;
        for iter_proto = 1:10
            th = zeros(1,50);
            B_record = zeros(50,(n1+i_ext)*n_v);
            parfor i = 1:50
                new_row = gen_random_row(n_v,M,extension_candidate,extension_prob);
                %%%% control error floor %%%%
                zero_pos = (new_row == 0) .* (B_bats_opt_old(end,:) == 0) .* (B_bats_opt_old(end-1,:) == 0);
                if (sum(zero_pos) > 0)
                    zero_id = find(zero_pos == 1);               
                    new_row(zero_id(1)) = randi([1 2],1,1);
                end
                %%%% control error floor end %%%%
                B_bats = [B_bats_opt_old;new_row];
                B_record(i,:) = reshape(B_bats,1,(n1+i_ext)*n_v);
                th(i) = Proto_BATS_Threshold_Approximate(B_bats,B_p,M,q,zeta,punc_vector_opt,n_hops);
            end
            [th_max,th_max_id] = max(th);
            if (th_max > th_opt)
                find_a_new_protograph = 1;
                th_opt = th_max;
                B_bats_opt = reshape(B_record(th_max_id,:),n1+i_ext,n_v);                  
            end
        end
        if (find_a_new_protograph == 1)
            rate = (size(B_bats_opt,2)-size(B_p,1))/(size(B_bats_opt,1)-sum(punc_vector_opt));
            th_capacity = capacity(n_hops,th_opt*ones(1,n_hops),M,q,zeta);
            fprintf('Extension Rows = %d, iter = %d, th_erasure = %.4f, th_capacity = %.4f, Rate = %.4f\n',i_ext,iter,th_opt,th_capacity,rate);
        end
    end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate a random protogaph of size m*n with some constraints
% this function will not fail to generate a row of a given degree
function [B] = gen_random_protograph(m,n,row_deg,light_candidate,dark_candidate,light_prob,dark_prob,dark_area)
    B = zeros(m,n);
    for i = 1:m
        for j = 1:n-1
            remaining_max = sum(dark_area(i,j+1:end)*dark_candidate(end) + (1-dark_area(i,j+1:end))*light_candidate(end));
            current_min = max(row_deg(i) - sum(B(i,:)) - remaining_max,0);
            if (dark_area(i,j) == 1)
                candidate = dark_candidate(current_min+1:end);
                prob = dark_prob(current_min+1:end)/sum(dark_prob(current_min+1:end));
            else
                candidate = light_candidate(current_min+1:end);
                prob = light_prob(current_min+1:end)/sum(light_prob(current_min+1:end));
            end

            B(i,j) = randsrc(1,1,[candidate;prob]);

            if (sum(B(i,:)) > row_deg(i))
                B(i,j) = B(i,j) - (sum(B(i,:)) - row_deg(i));   % tends to generate more zeros at the end
            end          
            if (B(i,j) < 0)
                error('B(i,j) < 0!');
            end
        end
        B(i,n) = row_deg(i) - sum(B(i,:));
        % The following operations remove the tendency of zeros at the end
        light_vec = B(i,dark_area(i,:) == 0);
        dark_vec = B(i,dark_area(i,:) == 1);
        light_vec = light_vec(randperm(size(light_vec,2)));
        dark_vec = dark_vec(randperm(size(dark_vec,2)));
        B(i,dark_area(i,:) == 0) = light_vec;      
        B(i,dark_area(i,:) == 1) = dark_vec;
    end
end

% generate a row of degree <= dg_max
function [B] = gen_random_row(n,dg_max,candidate,prob)

B = zeros(1,n);
for i = 1:n
    B(i) = randsrc(1,1,[candidate;prob]);
    if (sum(B) > dg_max)
        B(i) = B(i) - (sum(B) - dg_max);
    end
end
B = B(randperm(n));

end

% When calculating the BATS-check to variable message, the error probability is averaged by all incoming messages rather than the exact extrinsic messages
% Here, due to Homogeneous Links, the threshold measured by "capacity" in the paper is equivalent to the threshold measured by "erasure probability"
function [th] = Proto_BATS_Threshold_Approximate(B_bats,B_p,M,q,zeta,punc_vector,network_length)

[n1,k] = size(B_bats);
dg_col = sum(B_bats,1);
dg_row = sum(B_bats,2);
[n2,~] = size(B_p);
n = n1+n2;
B = [B_bats;B_p];

e_up = 1;
e_low = 0;
while (e_up - e_low > 0.001)
e = (e_up + e_low)/2;
h = Get_LinearNet_h(network_length,e*ones(1,network_length),M,q,zeta);
hbar = get_hbar(M,h,q);

x = ones(n,k);
y = ones(n,k);
z = ones(1,k);
for it = 2:1000
    % variable-to-check update
    for j = 1:k
        tmp_prod = prod(y(:,j).^B(:,j));
        for i = 1:n
            if (B(i,j) <= 0)
                x(i,j) = 1;
            else
                if (y(i,j) <= 0)
                    tmp = 1;
                    for i_prime = 1:n
                        if (i_prime ~= i)
                            tmp = tmp * y(i_prime,j)^B(i_prime,j);
                        end
                    end
                    tmp = tmp * y(i,j)^(B(i,j)-1);
                    x(i,j) = tmp;
                else
                    x(i,j) = tmp_prod/y(i,j);
                end
            end
        end
    end

    % check-to-variable update
    for i = 1:n1
        xbar = (sum(B(i,:).*x(i,:)))/(dg_row(i));     
        y(i,:) = Psi(xbar,dg_row(i),hbar,M);
        y(i,:) = double(B(i,:)>0).*((1-punc_vector(i))*y(i,:)+punc_vector(i)) + double(B(i,:)==0);
    end
    for i = n1 + 1:n
        tmp_prod = prod((1-x(i,:)).^B(i,:));
        for j = 1:k
            if (B(i,j)<=0)
                y(i,j) = 1;
            else
                if (1-x(i,j) <= 0)
                    tmp = 1;
                    for j_prime = 1:k                    
                        if (j_prime ~= j)
                            tmp = tmp*((1-x(i,j_prime))^B(i,j_prime));
                        else
                            tmp = tmp*((1-x(i,j))^(B(i,j)-1));
                        end
                    end
                    y(i,j) = 1-tmp;
                else
                    y(i,j) = 1 - tmp_prod/(1-x(i,j));
                end
            end
        end
    end
    
    % APP update
    for j = 1:k
        z(j) = 1;
        for i = 1:n
            z(j) = z(j)*y(i,j)^B(i,j);
        end
    end

    % check convergence
    Converge = 1;
    for j = 1:k
        if (z(j) > 1e-4)
            Converge = 0;
            break;
        end
    end
    if (Converge == 1)
        break;
    end
end
    if (Converge == 1)
        e_low = e;
    else
        e_up = e;
    end
end

th = e_low;

end

function [y] = Psi(x,d,hbar,M)
    y = 0;
    for r = 1:min(d-1,M)
        y = y + betainc(1-x,d-r,r)*hbar(r+1);
    end
    for r = d:M
        y = y + hbar(r+1);
    end
    y = 1-y;
end

function [hbar] = get_hbar(M,h,q)
    hbar = zeros(1,M+1);
    for r = 0:M
        for s = r:M
            hbar(r+1) = hbar(r+1) + zeta_func(s,r,q)/q^(s-r)*h(s+1);
        end
    end
end

function [zeta_table] = generate_zeta(M,q)
    for i = 0:M
        for j = 0:M
            zeta_table(i+1,j+1) = zeta_func(i,j,q);
        end
    end
end

function [C] = capacity(length,e,M,q,zeta)
h = Get_LinearNet_h(length,e,M,q,zeta);
C = sum((0:M).*h);
end
function [num_threads,BATS_type,num_input_sym,num_inter_sym,M,q,N_INAC,Psi,precode_G_sys,precode_H,BATS_Profile,k,e,n_array] = read_profile(file_name)
    Profile = fopen(file_name,'r');
    fscanf(Profile, '%s', 2);
    fscanf(Profile, '%s', 1);
    num_threads = fscanf(Profile, '%d', 1);
    fscanf(Profile, '%s', 1);
    BATS_type = fscanf(Profile, '%s', 1);
    fscanf(Profile, '%s', 1);
    num_input_sym = fscanf(Profile, '%d', 1);
    fscanf(Profile, '%s', 1);
    num_inter_sym = fscanf(Profile, '%d', 1);
    fscanf(Profile, '%s', 1);
    M = fscanf(Profile, '%d', 1);
    fscanf(Profile, '%s', 1);
    q = fscanf(Profile, '%d', 1);
    fscanf(Profile, '%s', 1);
    N_INAC = fscanf(Profile, '%d', 1);
    fscanf(Profile, '%s', 1);
    num_effective_dg = fscanf(Profile, '%d', 1);
    fscanf(Profile, '%s', 1);
    eff_dg = fscanf(Profile, '%d', num_effective_dg);
    eff_dg = eff_dg';   % change to row vector
    fscanf(Profile, '%s', 1);
    eff_dg_ratio = fscanf(Profile, '%f', num_effective_dg);
    eff_dg_ratio = eff_dg_ratio';
    fscanf(Profile, '%s', 1);
    precode_profile = fscanf(Profile, '%s', 1);
    fscanf(Profile, '%s', 1);
    precode_rank_option = fscanf(Profile, '%d', 1);
    fscanf(Profile, '%s', 1);
    batch_profile = fscanf(Profile, '%s', 1);
    fscanf(Profile, '%s', 1);
    k = fscanf(Profile, '%d', 1);
    fscanf(Profile, '%s', 1);
    e = fscanf(Profile, '%f', k);
    e = e';
    fscanf(Profile, '%s', 1);
    n_min = fscanf(Profile, '%d', 1);
    fscanf(Profile, '%s', 1);
    n_step = fscanf(Profile, '%d', 1);
    fscanf(Profile, '%s', 1);
    n_max = fscanf(Profile, '%d', 1);
    n_array = n_min:n_step:n_max;
    fclose(Profile);

    % initialization
    Psi = zeros(1,num_inter_sym);
    Psi(eff_dg) = eff_dg_ratio;
    Psi = Psi/sum(Psi);
    
    m = num_inter_sym - num_input_sym;
    if (num_input_sym > num_inter_sym)
        error('the number of intermediate symbols needs to >= the number of input symbols')
    elseif (num_input_sym == num_inter_sym)
        H_gf = gf([],log2(q));   % empty matrix        
    else
        LDPC_H = gf(cell2mat(struct2cell(load(precode_profile))),log2(q));
        if (precode_rank_option == 0)
            H_gf = LDPC_H;
            for i = 1:num_inter_sym-num_input_sym
                for j = 1:num_inter_sym
                    if (H_gf(i,j)~=0)
                        H_gf(i,j) = gf(1,log2(q));
                    end
                end
            end
        elseif (precode_rank_option == 1)
            H_gf = LDPC_H;
            for i = 1:num_inter_sym-num_input_sym
                for j = 1:num_inter_sym
                    if (H_gf(i,j)~=0)
                        H_gf(i,j) = gf(randi([1 q-1],1,1),log2(q));
                    end
                end
            end
        else
            while (1)
                H_gf = LDPC_H;
                for i = 1:num_inter_sym-num_input_sym
                    for j = 1:num_inter_sym
                        if (H_gf(i,j)~=0)
                            H_gf(i,j) = gf(randi([1 q-1],1,1),log2(q));
                        end
                    end
                end
                if (rank(H_gf)==min(size(H_gf)))
                    break;
                end
            end
        end
    end
    [precode_H,~,precode_G_sys,pos] = Generate_LDPC_G_Matrix_From_H_Matrix(H_gf,q,m,num_inter_sym);
    
    if (strcmp(BATS_type,'protograph'))
        batch_file = fopen(batch_profile,'r');
        num_batch = fscanf(batch_file,'%d',1);
        batch_matrix = zeros(num_batch,num_inter_sym);
        for i = 1:num_batch
            dg = fscanf(batch_file,'%d',1);        % degree
            tmp = fscanf(batch_file,'%d',dg);
            tmp = tmp';
            batch_matrix(i,tmp) = 1;
        end
        fclose(batch_file);
    
        batch_matrix = batch_matrix(:,pos);      % adjust the batch connections because the precode matrix's columns may be adjusted during Gaussian elimination
        BATS_Profile = cell(1,num_batch);
        for i = 1:num_batch
            BATS_Profile{1,i} = find(batch_matrix(i,:)>0);
        end
    else
        BATS_Profile = [];
    end
end
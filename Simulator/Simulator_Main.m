clear;clc;

[num_threads,BATS_Type,K_pure,K,M,q,N_INAC,Psi,precode_G_sys,precode_H,...
    BATS_Profile,k,e,n_array] = read_profile('Profile.txt');

delete(gcp('nocreate'));
parallel_obj = parpool(num_threads);
if (N_INAC == 0)
    Decoding_Method = 2;        % use fast BP decoder                                      
else
    Decoding_Method = 1;        % use real BP-INA decoder
end
m = K - K_pure;                 % the number of precoded batches

% Some initial variables
B = gf(zeros(1,K),log2(q));
FER = zeros(size(n_array));
SER = zeros(size(n_array));
inac_num_max = zeros(size(n_array));
inac_num_avg = zeros(size(n_array));
error_num_parallel = zeros(1,10*num_threads);
error_sym_num_parallel = zeros(1,10*num_threads);
inac_num_parallel = zeros(1,10*num_threads);


for n_index = 1:max(size(n_array))
    test_num = 0;
    error_num = 0;
    error_sym_num = 0;
    n = n_array(n_index);
    while (error_num < 100 || test_num < 100)
        parfor ii = 1:10*num_threads
        %for ii = 1:10*num_threads   
            %B_pure = gf(randi([0 q-1],T,K_pure),log2(q));
            %B = B_pure*precode_G_sys;
            B = gf(zeros(1,K),log2(q));     % transmit zero codeword

            [A,G,H,Y] = BATS_Transmission(B,M,q,Psi,k,e,N_INAC,precode_H,1,n,m,BATS_Type,BATS_Profile);

            if (Decoding_Method == 1)
                [B_dec,fail,inac_num] = BP_Decoder(A,G,H,Y,q,K,n,m,M,N_INAC);
                error_num_parallel(ii) = fail;
                error_sym_num_parallel(ii) = sum(sum(B~=B_dec));
                inac_num_parallel(ii) = inac_num;      
            else
                [num_undecoded] = fast_BP_decoder(A,G,H,q,K,n,m);
                error_num_parallel(ii) = double(num_undecoded > 0);
                error_sym_num_parallel(ii) = num_undecoded;
                inac_num_parallel(ii) = 0;
            end
        end
        if (max(inac_num_parallel) > inac_num_max(n_index))
            inac_num_max(n_index) = max(inac_num_parallel);
        end
        inac_num_avg(n_index) = inac_num_avg(n_index) + sum(inac_num_parallel);
        error_num = error_num + sum(error_num_parallel);
        error_sym_num = error_sym_num + sum(error_sym_num_parallel);
        test_num = test_num + 10*num_threads;
    end
    FER(n_index) = error_num/test_num
    SER(n_index) = error_sym_num/test_num/K
    inac_num_avg(n_index) = inac_num_avg(n_index)/test_num;
end


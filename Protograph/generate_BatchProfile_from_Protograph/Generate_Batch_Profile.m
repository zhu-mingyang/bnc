clear;clc;

% protograph of the BATS code B^(2) (excluding the protograph of precode)
B_bats = [
2	5	2	1	2	1	0	3
2	3	1	1	3	1	5	5
0	4	3	3	2	2	4	4
3	2	3	0	1	0	1	4
4	5	3	3	4	5	4	3
5	3	3	5	3	5	3	5
0	3	2	0	1	0	3	3
0	3	2	0	3	0	1	3
1	3	3	0	1	0	0	3
0	3	3	1	1	0	0	2
0	2	0	0	3	1	0	3
1	3	1	0	0	0	1	3
0	3	3	1	0	0	0	1
0	2	2	0	2	0	0	2
];
punc_vector = [0.870550872847841	0.930671513602045	0.926842168337794	0.930671513602045	0.960631965805137	0.960631965805137	0.940000000000000	0.940000000000000	0.940000000000000	0.940000000000000	0.940000000000000	0.940000000000000	0.940000000000000	0.940000000000000];

precode_H = cell2mat(struct2cell(load('./LDPC matrix/QC_LDPC_precode_960x2560_Z64.mat')));

Z1 = 5; Z2 = 64;      % two-step lifting factors
row_core = 6;         % the number of rows in the core protogmatrix (B_c^(2))
max_unconnected = 1000000;  % max number of unconnected VNs in the CORE protogmatrix. Here we set it to a very large value to disable it
min_unconnected = 650;      % min number of unconnected VNs in the CORE protogmatrix. This value may affect the performance. For the precode used in our paper (rate-5/8), this value = (0.2~0.3)*(the number of intermediate symbols) is good, but this is not a definitive conclusion.
n_rows_check_unconnection = 220;     % number of rows (upper submatrix) for unconnection ratio check

unconnected_ratio = generate_batch_profile(B_bats,row_core,punc_vector,Z1,Z2,precode_H,max_unconnected,min_unconnected,n_rows_check_unconnection);

% generate BATS connections
function [unconnected_ratio] = generate_batch_profile(B_bats,row_core,punc_vector,Z1,Z2,precode_H,max_punc,min_punc,n_rows_check_unconnection)

% get parameters of precode H matrix
[m,n] = size(precode_H);
dc = sum(precode_H,2);
edge = (precode_H > 0);

% lifting
[B_lift,~] = Lifting_Protograph(B_bats,Z1,Z2);

while (1)
    row_number_tot = 0;
    H = [];
    for i = 1:size(B_bats,1)
        row_number = ceil((1-punc_vector(i))*Z1*Z2);
        if (row_number <= 0)
            continue;
        end
        row_index = randperm(Z1*Z2,row_number);   % randomly choose (1-p_i)*Z1*Z2 rows
        %row_index = 1:row_number;                  % choose the first (1-p_i)*Z1*Z2 rows
        row_index = row_index + (i-1)*Z1*Z2;
        H(row_number_tot+1:row_number_tot+row_number,:) = B_lift(row_index,:);
        if (i == row_core)
            H_core = H;
        end
        row_number_tot = row_number_tot + row_number;
    end

    % check connection ratio of the upper submatrix of B^(2) (We hope that some symbols remain disconnected by B^(2), as they can be resolved by the precode.)
    d = zeros(1,n);
    d(sum(H(1:n_rows_check_unconnection,:))>0) = 1;
    if (sum(d==0)>max_punc || sum(d==0)<min_punc)
        %fprintf('unconnection ratio unsatisfied, continuous to generate BATS codes...\n');
        continue;
    end
    unconnected_ratio = sum(d==0)/n;  

    % check BP decoding
    d = zeros(1,n);
    d(sum(H_core)>0) = 1;
    while (1)
        new_dec_vn = 0;
        for i = 1:m
            if (sum(d(edge(i,:))) == dc(i)-1)
                d(edge(i,:)) = 1;
                new_dec_vn = new_dec_vn + 1;
            end
        end
        if (new_dec_vn == 0)
            break;
        end
    end
    if (sum(d==0)==0)   % precode decoding succeeds
        break;
    else
        fprintf('BP decoding fails, continuous to generate BATS codes...\n');
    end
end

BATS_Profile = cell(1,row_number_tot);
for i = 1:row_number_tot
    BATS_Profile{1,i} = find(H(i,:)>0);
end

% write to file
outfile = fopen('./batch profile/batch_profile.txt','w');
fprintf(outfile,'%d\n',row_number_tot);
for i = 1:row_number_tot
    fprintf(outfile,'%d ',length(BATS_Profile{1,i}));
    for j = 1:length(BATS_Profile{1,i})
        fprintf(outfile,'%d ',BATS_Profile{1,i}(j));
    end
    fprintf(outfile,'\n');
end
fclose(outfile);

end
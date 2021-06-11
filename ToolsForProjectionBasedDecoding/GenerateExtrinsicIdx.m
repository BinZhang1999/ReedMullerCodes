 function projection_extrinsic_idx = GenerateExtrinsicIdx(projection_idx)
 %% Generate the extrinsic idx for projection aggregation decoding
 %  Input:
 %      projection_idx: size(projection_idx) = [num_overlap, num_coset,
 %      num_subspace]; 
 %      !!! The zero idx in projection_idx is represented by 1 !!!
 %  Output:
 %      projection_extrinsic_idx: [num_overlap-1, N, num_subspace] =
 %      size(projection_extrinsic_idx);
 %%
       [num_overlap, num_coset, num_subspace] = size(projection_idx);
       N = num_coset * num_overlap;
       projection_extrinsic_idx = zeros(num_overlap-1, N, num_subspace);
       for i = 1:num_subspace
            for k = 1:num_coset
                for j = 1:num_overlap
                    projection_extrinsic_idx(1:end,projection_idx(j,k,i),i) = ...
                    setdiff(projection_idx(:,k,i),projection_idx(j,k,i));
                end
            end
       end 
 end
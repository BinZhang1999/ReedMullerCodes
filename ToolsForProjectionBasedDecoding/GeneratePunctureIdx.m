function puncture_idx = GeneratePunctureIdx(projection_idx)
%% Generate puncture_idx from projection_idx
% Input:
%   projection_idx: The projection_idx of CPA for RM(m-r-1,m).
%   [num_overlap, num_coset, num_subspace] = size(projection_idx);
% Output:
%   puncture_idx: [num_cn, num_coset] = size(puncture_idx)
%%
[num_overlap, num_coset, num_subspace] = size(projection_idx);
num_cn = num_overlap * num_subspace;
puncture_idx = zeros(num_cn , num_coset);

for i = 1:num_subspace
   for j = 1:num_overlap
      idx = (i-1)*num_overlap+j;
      puncture_idx(idx, :) = projection_idx(j,:,i); 
   end
end

end
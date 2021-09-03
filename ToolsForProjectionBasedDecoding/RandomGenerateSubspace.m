function subspace_list = RandomGenerateSubspace(m,k,num_subspace)
%% Randomly Generate subspace
% Input: 
%   m: m-dimensional EG over GF(2)
%   k: k-dimensional subspace
%   num_subspace: required subspace number. It should bigger than
%   nchoosek(n,k)
% Output:
%   subspace_list: [num_subspace, num_point] = size(subspace_list)
%                  num_subspace is the subspace number in subspace_list
%                  num_point is the point number in one subspace
%                  the points are represented by the dec form 
%                  eg: [(000),(001),(010),(011)] are represented by
%                      [0,1,2,3]
%                  the number is in a ascend manner from left to right
%%
% Init 
num_point = 2^k;
p_s = 0;
subspace_list = zeros(num_subspace , num_point);
% first generate the trivial subspace whose base vector has only one nonzero entry
base_vec_list = myDec2Bin([0:2^m-1]',m);
base_vec_list = base_vec_list(2:end,:);
base_idx = 2.^(0:m-1);
idx_list = nchoosek(base_idx, k);
u_list = myDec2Bin([0:2^k-1]',k);
for idx = idx_list'
    base_vec = base_vec_list(idx,:);
    subspace = mod(u_list * base_vec, 2);
    subspace_dec = sort(myBin2Dec(subspace, m),'ascend');
    p_s = p_s+1;
    subspace_list(p_s,:) = subspace_dec;
end

while p_s < num_subspace
    base_vec = randn(k,m)>0;
    if rank(gf(base_vec))<k % need to be full rank
        continue;
    end
    subspace = mod(u_list * base_vec, 2);
    subspace_dec = sort(myBin2Dec(subspace, m),'ascend');
    subspace_dec = subspace_dec';
    if any(ismember(subspace_list(1:p_s,:),subspace_dec,'rows'))
        continue; % need to be unique
    end
    
    p_s = p_s+1;
    subspace_list(p_s,:) = subspace_dec;
end

end
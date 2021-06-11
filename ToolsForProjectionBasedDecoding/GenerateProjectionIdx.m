function projection_idx = GenerateProjectionIdx(subspace_list, m, k)
%% Generate the projection_idx from the subspace list
% Input:
%   m,k: the subspace is k-dim of EG(m,2)
%   subspace_list: [num_subspace, num_point] = size(subspace_list)
%                  num_subspace is the subspace number in subspace_list
%                  num_point is the point number in one subspace
%                  the points are represented by the dec form 
%                  eg: [(000),(001),(010),(011)] are represented by
%                      [0,1,2,3]
%                  the number is in a ascend manner from left to right
% Ouput:
%   projection_idx: [num_overlap, num_coset, num_subspace] =
%                   size(projection_idx). 
%                   num_overlap is the number of the points which are
%                   projected onto the same position. 
%                   num_coset is the number of the cosets of one subspace. 
%                   It is also the length of the projected vector.
%                   num_subspace is the number of the subspaces.
%                   !!! the zero is represented by 1, so k is represented
%                   by k+1 !!!
% Eg:
% % Generate the projection idx from subspace_list_m7k2
% projection_idx_m7k2 = GenerateProjectionIdx(subspace_list_m7k2,7,2);
%%
% Init
N = 2^m; % the point num of EG(m,2)
all_points = 0:(N-1);

[num_subspace, num_overlap] = size(subspace_list);
num_coset = 2^(m-k);
projection_idx = zeros(num_overlap, num_coset, num_subspace);

% enumerate the cosets of each subspace
for i_subspace = 1:num_subspace
   p_coset = 1;
   projection_idx(:, p_coset, i_subspace) = subspace_list(i_subspace,:);
   points_uncontained = setdiff(all_points, subspace_list(i_subspace,:));
   
   subspace_bin = MyDec2Bin(subspace_list(i_subspace,:)',m);
   % caculate the coset
   for point = points_uncontained
      if any(ismember(projection_idx(:,:,i_subspace), point),'all')
          continue;
      end
      point_bin = MyDec2Bin(point, m);
      coset_bin = (subspace_bin~=point_bin);
      coset = MyBin2Dec(coset_bin, m);
      p_coset = p_coset+1;
      projection_idx(:,p_coset,i_subspace) = coset;
   end
end
projection_idx = projection_idx+1;
end
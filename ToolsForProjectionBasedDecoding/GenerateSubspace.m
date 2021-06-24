function subspace_list = GenerateSubspace(m, k, subspace_list_last)
%% Generate all the k-dim subspaces of EG(m,2) from (k-1) dim subspaces
% Input:
%   m: m-dimensional EG over GF(2)
%   k: k-dimensional subspace
%   subspace_list_last: [num_subspace_last, num_point_last] = size(subspace_list_last)
%                  num_subspace_last is the subspace number in
%                  subspace_list_last
%                  num_point_last is the point number in one subspace
%                  the points are represented by the dec form 
%                  eg: [(000),(001),(010),(011)] are represented by
%                      [0,1,2,3]
%                  the number is in a ascend manner from left to right
% Ouput: 
%   subspace_list: [num_subspace, num_point] = size(subspace_list)
%                  num_subspace is the subspace number in subspace_list
%                  num_point is the point number in one subspace
%                  the points are represented by the dec form 
%                  eg: [(000),(001),(010),(011)] are represented by
%                      [0,1,2,3]
%                  the number is in a ascend manner from left to right
% Eg:
% % Generate all the 1-dimensional subspaces if EG(7,2) 
% subspace_list_m7k1 = GenerateSubspace(7, 1, [0]);
% % Generate all the 2-dimensional subspaces if EG(7,2) 
% subspace_list_m7k2 = GenerateSubspace(7, 2, subspace_list_m7k1);

%%
% Init
N = 2^m; % the point num of EG(m,2)
all_points = 0:(N-1);

[num_subspace_last, ~] = size(subspace_list_last);

num_subspace = CountSubspace(m,k,1);
num_point = 2^k;
% sometimes need to replace the num_subsapce and num_point to the real value
% eg: subspace_list = zeros(11811, 16);
subspace_list = zeros(num_subspace, num_point);
p_subspace = 0; % the pointer

% add the new point to the last subspace
for i_subspace_last = 1:num_subspace_last
% for i_subspace_last = num_subspace_last:-1:1
    subspace_last = subspace_list_last(i_subspace_last,:);
    points_uncontained = setdiff(all_points, subspace_last);
    subspace_last_bin = MyDec2Bin(subspace_last', m);
    
    for point = points_uncontained
        % add a new point to the (k-1)-dim subspace and generate a k-dim
        % subspace
        point_bin = MyDec2Bin(point, m);
        points_add_bin = (subspace_last_bin ~= point_bin);
        points_add = MyBin2Dec(points_add_bin, m);    
        subspace = sort([subspace_last, points_add'],'ascend');
        
        % if the tmp_subspace is already in the subspace list
        if any(ismember(subspace_list(1:p_subspace,:), subspace, 'rows'))
            continue; 
        end
        
        % if not, add the new subspace to the list
        p_subspace = p_subspace+1;
        subspace_list(p_subspace,:) = subspace;
        if p_subspace == num_subspace
            return;
        end
    end
end


end
function idx_set = IndexByInfoBit(projection_idx, G, r, m)
%% Index the projection_idx by info bits 
% Input: 
%   projection_idx: [num_overlap, num_coset,
%   num_subspace]=size(projection_idx);
%   G: The generator matrix of RM code. 
%   r,m: RM(r,m) code
% Ouput:
%   idx_set: [K_rth, max_num_subspace_infobit+1] = size(idx_set)
%            the first number is the number of the num_subspace of this    
%            info bit
%%

G_rth = G(sum(G,2)==2^(m-r),:);
[num_overlap, num_coset, num_subspace]= size(projection_idx);
K_rth = nchoosek(m,r);
idx_set = zeros(K_rth, 2^(r*(m-r))+1);
cnt = zeros(K_rth,1);

for i = 1:num_subspace
   G_x = zeros(K_rth, num_coset); 
   for j = 1:num_overlap
      G_x = G_x + G_rth(:, projection_idx(j,:,i)); 
   end
   G_x = mod(G_x,2);
   
   for k = 1:K_rth
      if G_x(k,1)
          cnt(k) = cnt(k)+1;
          idx_set(k,cnt(k)) = i;
      end
   end
end

for k = 1:K_rth
    idx_nonzero = (idx_set(k,:)~=0);
    num_subspace_this_infobit = sum(idx_nonzero);
    idx_set(k,2:(num_subspace_this_infobit+1)) = sort(idx_set(k,idx_nonzero),'ascend');
    idx_set(k,1) = num_subspace_this_infobit;
end

max_num_subspace_infobit = max(idx_set(:,1));
idx_set = idx_set(:,1:(max_num_subspace_infobit+1));

end
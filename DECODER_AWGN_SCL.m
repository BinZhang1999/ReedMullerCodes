classdef DECODER_AWGN_SCL
properties
    L; % list size 
    K; % info bits num
    N;
    m; % N = 2^m
    frozen_bits; % the vector 1-frozen bit 0-nonfrozen bit
    lambda_offset;
    llr_layer_vec;
    bit_layer_vec;
    
    Gm;
    G;
end

methods
    function obj = Init(obj, L, K, N, frozen_bits)
       obj.L = L;
       obj.K = K;
       obj.m = log2(N);
       obj.N = N;
       obj.frozen_bits = frozen_bits;
       obj.lambda_offset = 2.^(0 : log2(N));
       obj.llr_layer_vec = get_llr_layer(N);
       obj.bit_layer_vec = get_bit_layer(N);
       
       F = [[1 0]
            [1 1]];
       obj.Gm = F;
       for i = 2:obj.m
           obj.Gm = kron(F, obj.Gm);
       end
       obj.G = obj.Gm(~obj.frozen_bits,:);
     
    end
    
    function [u_hat, v_hat] = Decode(obj, rx, sigma)
        u_hat = zeros(1,obj.K)-1;
        v_hat = zeros(1,obj.N)-1;
        llr = 2*rx/(sigma^2);
        
        lazy_copy = zeros(obj.m, obj.L);%If Lazy Copy is used, there is no data-copy in the decoding process. We only need to record where the data come from. Here,data refer to LLRs and partial sums.
        %Lazy Copy is a relatively sophisticated operation for new learners of polar codes. If you do not understand such operation, you can directly copy data.
        %If you can understand lazy copy and you just start learning polar codes
        %for just fews days, you are very clever,
        P = zeros(obj.N - 1, obj.L); %Channel llr is public-used, so N - 1 is enough. 用来存储信道译码过程中的llr 我们用最后一个判决
        C = zeros(obj.N - 1, 2 * obj.L);%I do not esitimate (x1, x2, ... , xN), so N - 1 is enough.
        u = zeros(obj.K, obj.L);%unfrozen bits that polar codes carry, including crc bits.
        PM = zeros(obj.L, 1);%Path metrics
        activepath = zeros(obj.L, 1);%Indicate if the path is active. '1'→active; '0' otherwise.
        cnt_u = 1;%information bit counter 
        %initialize
        activepath(1) = 1;
        lazy_copy(:, 1) = 1;
        %decoding starts
        %default: in the case of path clone, the origianl path always corresponds to bit 0, while the new path bit 1.
        for phi = 0 : obj.N - 1
            layer = obj.llr_layer_vec(phi + 1);
            phi_mod_2 = mod(phi, 2);
            for l_index = 1 : obj.L
                if activepath(l_index) == 0 % 如果路径没有激活  直接跳过
                    continue;
                end
                switch phi%Decoding bits u_0 and u_N/2 needs channel LLR, so the decoding of them is separated from other bits. 
                case 0
                index_1 = obj.lambda_offset(obj.m); % 正在第 m 层
                for beta = 0 : index_1 - 1
%                     zb = sign(llr(beta + 1)) * sign(llr(beta + index_1 + 1)) * min(abs(llr(beta + 1)), abs(llr(beta + index_1 + 1)));
                    P(beta + index_1, l_index) = log( (1+exp(llr(beta + 1)+llr(beta + index_1 + 1))) / (exp(llr(beta + 1)) + exp(llr(beta + index_1 + 1))) );
                end
                for i_layer = obj.m - 2 : -1 : 0   % 从 m-1 层开始往后算
                    index_1 = obj.lambda_offset(i_layer + 1);  % 这一层的 llr 的数量
                    index_2 = obj.lambda_offset(i_layer + 2);  % 上一层的 llr 的数量
                    for beta = 0 : index_1 - 1
%                         zb = sign(P(beta + index_2, l_index)) *...
%                             sign(P(beta + index_1 + index_2, l_index)) * min(abs(P(beta + index_2, l_index)), abs(P(beta + index_1 + index_2, l_index)));
                        P(beta + index_1, l_index) = log( (1+exp(P(beta + index_2, l_index)+P(beta + index_1 + index_2, l_index))) / (exp(P(beta + index_2, l_index)) + exp(P(beta + index_1 + index_2, l_index))) );
                    end
                end
                case obj.N/2
                index_1 = obj.lambda_offset(obj.m);
                for beta = 0 : index_1 - 1
                    x_tmp = C(beta + index_1, 2 * l_index - 1);
                    P(beta + index_1, l_index) = (1 - 2 * x_tmp) * llr(beta + 1) + llr(beta + 1 + index_1);
                end
                for i_layer = obj.m - 2 : -1 : 0
                    index_1 = obj.lambda_offset(i_layer + 1);
                    index_2 = obj.lambda_offset(i_layer + 2);
                    for beta = 0 : index_1 - 1
%                         zb = sign(P(beta + index_2, l_index)) *...
%                             sign(P(beta + index_1 + index_2, l_index)) * min(abs(P(beta + index_2, l_index)), abs(P(beta + index_1 + index_2, l_index)));
                        P(beta + index_1, l_index) = log( (1+exp(P(beta + index_2, l_index)+P(beta + index_1 + index_2, l_index))) / (exp(P(beta + index_2, l_index)) + exp(P(beta + index_1 + index_2, l_index))) );
                    end
                end
                otherwise
                index_1 = obj.lambda_offset(layer + 1);
                index_2 = obj.lambda_offset(layer + 2);
                for beta = 0 : index_1 - 1
                    P(beta + index_1, l_index) = (1 - 2 * C(beta + index_1, 2 * l_index - 1)) * P(beta + index_2, lazy_copy(layer + 2, l_index)) +...
                        P(beta + index_1 + index_2, lazy_copy(layer + 2, l_index));
                end
                for i_layer = layer - 1 : -1 : 0 % 总是一次 g 运算后 跟 f  运算到底
                    index_1 = obj.lambda_offset(i_layer + 1);
                    index_2 = obj.lambda_offset(i_layer + 2);
                    for beta = 0 : index_1 - 1
%                         zb = sign(P(beta + index_2, l_index)) *...      
%                             sign(P(beta + index_1 + index_2, l_index)) * min(abs(P(beta + index_2, l_index)),...
%                             abs(P(beta + index_1 + index_2, l_index)));
                        P(beta + index_1, l_index) = log( (1+exp(P(beta + index_2, l_index)+P(beta + index_1 + index_2, l_index))) / (exp(P(beta + index_2, l_index)) + exp(P(beta + index_1 + index_2, l_index)) ));
                    end
                end
        end
    end
    if obj.frozen_bits(phi + 1) == 0%if now we decode an unfrozen bit
        PM_pair = realmax * ones(2, obj.L); % 路径度量修改值的初始化  直接干到 inf
        % 第 l 列就表示 第 l list的情况
        % 第 1 / 2 行分别对应译码为 0 / 1 的路径值
        for l_index = 1 : obj.L
            if activepath(l_index) == 0
                continue;
            end
            if P(1, l_index) >= 0   % 判定为 0 就+llr值
                PM_pair(1, l_index) = PM(l_index);
                PM_pair(2, l_index) = PM(l_index) + P(1, l_index);
            else  % 判定为 1 就-llr值
                PM_pair(1, l_index) = PM(l_index) - P(1, l_index);
                PM_pair(2, l_index) = PM(l_index);
            end
        end
        middle = min(2 * sum(activepath), obj.L); % 最终保留的路径的数目   
        PM_sort = sort(PM_pair(:)); % 从小到大 给路径的度量排序   此时middle 也就代表着保留路径的最后一个
        PM_cv = PM_sort(middle); % 保留下来的最大的路径
        compare = PM_pair <= PM_cv; % 找出PM_pair保留下来的路径 compare（i） = 1 就保留   此时注意  可能会有大于等于 middle 数量的路径被标志为 1 因为可能PM_sort(middle) 有重复的情况
        
        % 避免出现middle位置有重复值的情况
        % 若出现 此时只保留 middle 个compair 为 1 
        if sum(sum(compare))~=middle
            while sum(sum(compare))~=middle
                [tmp_idx,~] = sort(PM_pair(:));
                compare(tmp_idx(middle+1))=0;
            end
        end
        
        
        kill_index = zeros(obj.L, 1);%to record the index of the path that is killed
        kill_cnt = 0;%the total number of killed path
        %the above two variables consist of a stack
        for i = 1 : obj.L
            if (compare(1, i) == 0)&&(compare(2, i) == 0)%which indicates that this path should be killed
                % 这意味着这一列路径都不行  都会被删除
                activepath(i) = 0;% 将这个路径结束
                kill_cnt = kill_cnt + 1;%push stack  计算一下kill_cnt的数目  前面说过 compare 里面标志为 1 的会多  因此标志为 0 的会少 所以kill_cnt 会偏少 所以后续会存在 kill_cnt 等于 0的情况 在出栈的时候
                kill_index(kill_cnt) = i; 
            end
        end
        % 上面这个for循环  列举了所有路径：译码为 0 和 1 都会被砍掉的情况
        for l_index = 1 : obj.L
            if activepath(l_index) == 0
                continue;
            end
            path_state = compare(1, l_index) * 2 + compare(2, l_index); % 这条路径的状态  对应的是  极化码讲义里的三种状态
            switch path_state%path_state can equal to 0, but in this case we do no operation.
                case 1 % 只有译码为 1 的路径得到保留
                    u(cnt_u, l_index) = 1;
                    C(1, 2 * l_index - 1 + phi_mod_2) = 1;
                    PM(l_index) = PM_pair(2, l_index);
                case 2 % 只有译码为 0 的路径得到保留
                    u(cnt_u, l_index) = 0;
                    C(1, 2 * l_index - 1 + phi_mod_2) = 0;
                    PM(l_index) = PM_pair(1, l_index);
                case 3 % 译码为 0 和译码为 1 的路径都得到保留 此时需要占用之前kill掉的路径的位置
                    index = kill_index(kill_cnt); % 把kill cnt里面的第一个路径位置拿出来
                    kill_cnt = kill_cnt - 1;%pop stack  左边即为堆栈的顶部
                    activepath(index) = 1; % 激活这个被再次利用的 kill path
                    %lazy copy
                    lazy_copy(:, index) = lazy_copy(:, l_index); % 记录一下此时  这条路径的lazy copy 即它从哪儿来的（与当下路径共有此前的数据，此后分道扬镳）
                    u(:, index) = u(:, l_index);
                    u(cnt_u, l_index) = 0;
                    u(cnt_u, index) = 1;
                    C(1, 2 * l_index - 1 + phi_mod_2) = 0;
                    C(1, 2 * index - 1 + phi_mod_2) = 1;
                    PM(l_index) = PM_pair(1, l_index);
                    PM(index) = PM_pair(2, l_index);
            end
        end
        cnt_u = cnt_u + 1;
    else%frozen bit operation
        for l_index = 1 : obj.L
            if activepath(l_index) == 0 % 对未激活的路径无需操作
                continue;
            end
            if P(1, l_index) < 0
                PM(l_index) = PM(l_index) - P(1, l_index); % 如果判决比特与frozen bit既定值不同  那么算路径度量时  要减去这个位置的llr值
            end
            if phi_mod_2 == 0  % 对于编号偶数的叶子节点  需要更新中间的译码比特
                C(1, 2 * l_index - 1) = 0;
            else
                C(1, 2 * l_index) = 0;
            end 
        end
    end 
    
    for l_index = 1 : obj.L%partial-sum return % 去返回一些部分的中间的译码结果
        if activepath(l_index) == 0
            continue
        end
        if (phi_mod_2  == 1) && (phi ~= obj.N - 1)
            layer = obj.bit_layer_vec(phi + 1);
            for i_layer = 0 : layer - 1
                index_1 = obj.lambda_offset(i_layer + 1);
                index_2 = obj.lambda_offset(i_layer + 2);
                for beta = index_1 : 2 * index_1 - 1
                    C(beta + index_1, 2 * l_index) = mod(C(beta, 2 *  lazy_copy(i_layer + 1, l_index) - 1) + C(beta, 2 * l_index), 2);%Left Column lazy copy
                    C(beta + index_2, 2 * l_index) = C(beta, 2 * l_index);   
                end
            end
            index_1 = obj.lambda_offset(layer + 1);
            index_2 = obj.lambda_offset(layer + 2);
            for beta = index_1 : 2 * index_1 - 1
                C(beta + index_1, 2 * l_index - 1) = mod(C(beta, 2 * lazy_copy(layer + 1, l_index) - 1) + C(beta, 2 * l_index), 2);%Left Column lazy copy
                C(beta + index_2, 2 * l_index - 1) = C(beta, 2 * l_index);
            end 
        end
    end
    %lazy copy
    if phi < obj.N - 1
        for i_layer = 1 : obj.llr_layer_vec(phi + 2) + 1
            for l_index = 1 : obj.L
                lazy_copy(i_layer, l_index) = l_index;
            end
        end
    end
end
%path selection.
[~, path_ordered] = sort(PM);
% 路径的顺序
u_hat = u(:,path_ordered(1))';
v_hat = mod(u_hat*obj.G,2);
        
    end
end
end
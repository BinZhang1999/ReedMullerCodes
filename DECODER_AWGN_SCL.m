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
        P = zeros(obj.N - 1, obj.L); %Channel llr is public-used, so N - 1 is enough. �����洢�ŵ���������е�llr ���������һ���о�
        C = zeros(obj.N - 1, 2 * obj.L);%I do not esitimate (x1, x2, ... , xN), so N - 1 is enough.
        u = zeros(obj.K, obj.L);%unfrozen bits that polar codes carry, including crc bits.
        PM = zeros(obj.L, 1);%Path metrics
        activepath = zeros(obj.L, 1);%Indicate if the path is active. '1'��active; '0' otherwise.
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
                if activepath(l_index) == 0 % ���·��û�м���  ֱ������
                    continue;
                end
                switch phi%Decoding bits u_0 and u_N/2 needs channel LLR, so the decoding of them is separated from other bits. 
                case 0
                index_1 = obj.lambda_offset(obj.m); % ���ڵ� m ��
                for beta = 0 : index_1 - 1
%                     zb = sign(llr(beta + 1)) * sign(llr(beta + index_1 + 1)) * min(abs(llr(beta + 1)), abs(llr(beta + index_1 + 1)));
                    P(beta + index_1, l_index) = log( (1+exp(llr(beta + 1)+llr(beta + index_1 + 1))) / (exp(llr(beta + 1)) + exp(llr(beta + index_1 + 1))) );
                end
                for i_layer = obj.m - 2 : -1 : 0   % �� m-1 �㿪ʼ������
                    index_1 = obj.lambda_offset(i_layer + 1);  % ��һ��� llr ������
                    index_2 = obj.lambda_offset(i_layer + 2);  % ��һ��� llr ������
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
                for i_layer = layer - 1 : -1 : 0 % ����һ�� g ����� �� f  ���㵽��
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
        PM_pair = realmax * ones(2, obj.L); % ·�������޸�ֵ�ĳ�ʼ��  ֱ�Ӹɵ� inf
        % �� l �оͱ�ʾ �� l list�����
        % �� 1 / 2 �зֱ��Ӧ����Ϊ 0 / 1 ��·��ֵ
        for l_index = 1 : obj.L
            if activepath(l_index) == 0
                continue;
            end
            if P(1, l_index) >= 0   % �ж�Ϊ 0 ��+llrֵ
                PM_pair(1, l_index) = PM(l_index);
                PM_pair(2, l_index) = PM(l_index) + P(1, l_index);
            else  % �ж�Ϊ 1 ��-llrֵ
                PM_pair(1, l_index) = PM(l_index) - P(1, l_index);
                PM_pair(2, l_index) = PM(l_index);
            end
        end
        middle = min(2 * sum(activepath), obj.L); % ���ձ�����·������Ŀ   
        PM_sort = sort(PM_pair(:)); % ��С���� ��·���Ķ�������   ��ʱmiddle Ҳ�ʹ����ű���·�������һ��
        PM_cv = PM_sort(middle); % ��������������·��
        compare = PM_pair <= PM_cv; % �ҳ�PM_pair����������·�� compare��i�� = 1 �ͱ���   ��ʱע��  ���ܻ��д��ڵ��� middle ������·������־Ϊ 1 ��Ϊ����PM_sort(middle) ���ظ������
        
        % �������middleλ�����ظ�ֵ�����
        % ������ ��ʱֻ���� middle ��compair Ϊ 1 
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
                % ����ζ����һ��·��������  ���ᱻɾ��
                activepath(i) = 0;% �����·������
                kill_cnt = kill_cnt + 1;%push stack  ����һ��kill_cnt����Ŀ  ǰ��˵�� compare �����־Ϊ 1 �Ļ��  ��˱�־Ϊ 0 �Ļ��� ����kill_cnt ��ƫ�� ���Ժ�������� kill_cnt ���� 0����� �ڳ�ջ��ʱ��
                kill_index(kill_cnt) = i; 
            end
        end
        % �������forѭ��  �о�������·��������Ϊ 0 �� 1 ���ᱻ���������
        for l_index = 1 : obj.L
            if activepath(l_index) == 0
                continue;
            end
            path_state = compare(1, l_index) * 2 + compare(2, l_index); % ����·����״̬  ��Ӧ����  �����뽲���������״̬
            switch path_state%path_state can equal to 0, but in this case we do no operation.
                case 1 % ֻ������Ϊ 1 ��·���õ�����
                    u(cnt_u, l_index) = 1;
                    C(1, 2 * l_index - 1 + phi_mod_2) = 1;
                    PM(l_index) = PM_pair(2, l_index);
                case 2 % ֻ������Ϊ 0 ��·���õ�����
                    u(cnt_u, l_index) = 0;
                    C(1, 2 * l_index - 1 + phi_mod_2) = 0;
                    PM(l_index) = PM_pair(1, l_index);
                case 3 % ����Ϊ 0 ������Ϊ 1 ��·�����õ����� ��ʱ��Ҫռ��֮ǰkill����·����λ��
                    index = kill_index(kill_cnt); % ��kill cnt����ĵ�һ��·��λ���ó���
                    kill_cnt = kill_cnt - 1;%pop stack  ��߼�Ϊ��ջ�Ķ���
                    activepath(index) = 1; % ����������ٴ����õ� kill path
                    %lazy copy
                    lazy_copy(:, index) = lazy_copy(:, l_index); % ��¼һ�´�ʱ  ����·����lazy copy �������Ķ����ģ��뵱��·�����д�ǰ�����ݣ��˺�ֵ�����
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
            if activepath(l_index) == 0 % ��δ�����·���������
                continue;
            end
            if P(1, l_index) < 0
                PM(l_index) = PM(l_index) - P(1, l_index); % ����о�������frozen bit�ȶ�ֵ��ͬ  ��ô��·������ʱ  Ҫ��ȥ���λ�õ�llrֵ
            end
            if phi_mod_2 == 0  % ���ڱ��ż����Ҷ�ӽڵ�  ��Ҫ�����м���������
                C(1, 2 * l_index - 1) = 0;
            else
                C(1, 2 * l_index) = 0;
            end 
        end
    end 
    
    for l_index = 1 : obj.L%partial-sum return % ȥ����һЩ���ֵ��м��������
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
% ·����˳��
u_hat = u(:,path_ordered(1))';
v_hat = mod(u_hat*obj.G,2);
        
    end
end
end
function [funs] = student_sols()
    
%     function tre = polynomial2trellis(ConstraintLength,Generator)
%         % Generator is a decimal value array
%         num_register = ConstraintLength - 1;
%         [a,b]= size(Generator);
%         tre.numInputSymbols = 2.^a;
%         tre.numOutputSymbols = 2.^b;
%         tre.numStates = 2.^(num_register);
%         % Initiate nextStates and outputs matrix
%         tre.nextStates = zeros(tre.numStates,tre.numInputSymbols);
%         tre.outputs = zeros(tre.numStates,tre.numInputSymbols);
%         
%         for i = 1:tre.numStates
%             for j = 1:tre.numInputSymbols
%                 d = [(int2bit(j-1,a))' (int2bit(i-1,num_register))'];
%                 tre.nextStates(i,j) = bi2de(d(end-num_register:end-1),'left-msb');
%                 
%                 e = d.*(int2bit(Generator(1),num_register+1))';
%                 f = d.*(int2bit(Generator(2),num_register+1))';
%                 g = mod([sum(e) sum(f)],2);
%                 tre.outputs(i,j) = bi2de(g,'left-msb');
%   
%             end
% 
%         end
%         
%     end
    
    function tre = polynomial2trellis(varargin)
        % Generator is a decimal value array
        [~,column] = size(varargin{1});
        num_register = sum(varargin{1}) - column;
        [a,b]= size(varargin{2});
        tre.numInputSymbols = 2.^a;
        tre.numOutputSymbols = 2.^b;
        tre.numStates = 2.^(num_register);
        % Initiate nextStates and outputs matrix
        tre.nextStates = zeros(tre.numStates,tre.numInputSymbols);
        tre.outputs = zeros(tre.numStates,tre.numInputSymbols);
        
        % create a new array to store the index of inserting input
        insert_index = ones(1,a);  
        for k = 1:column-1
            insert_index(k+1) = sum(varargin{1}(end-k+1:end))-k+1;
        end

        % Trellis Structure for Feedforward Convolutional Encoder
        if nargin == 2
            % i is row, j is column, i-1 is states, j-1 is input
            for i = 1:tre.numStates
                    % the state before update, shift 1 bit right, transform
                    % to binary
                    current_state = int2bit(i-1,num_register);
                    d = int2bit(bitshift(i-1,-1),num_register);
                for j = 1:tre.numInputSymbols
                    % update nextStates value
                    input = int2bit(j-1,a);               
                    if a > 1
                        d(1) = input(end);
                        for k = 1:column-1
                            d(insert_index(k+1)) = input(end-k);
                        end
                    else 
                        d(1) = input(end);
                    end
                    tre.nextStates(i,j) = bi2de(d','left-msb');
                    
                    % update outputs value
                    
                    generator_matrix = zeros(sum(varargin{1}),b);
                    state_input_matrix = zeros(sum(varargin{1}),1);
                    for m = 1:a
                       row_matrix = int2bit(varargin{2}(m,:),varargin{1}(m));
                       %d1 = [input(m);int2bit(i-1,num_register)(:)]
                       %lower = sum(varargin{1}(1:m-1));
                       %test = varargin{1}(m);
                       %generator_matrix(end-lower-varargin{1}(m)+1:end-lower,:) = row_matrix;
                       if m == 1
                            generator_matrix(end-varargin{1}(m)+1:end,:) = row_matrix;
                            state_input_matrix(end-varargin{1}(m)+1:end) = [input(m);current_state(insert_index(end):end)];
                       else  
                            lower = sum(varargin{1}(1:m-1));
                            generator_matrix(end-lower-varargin{1}(m)+1:end-lower,:) = row_matrix;
                            state_input_matrix(end-lower-varargin{1}(m)+1:end-lower) = [input(m);current_state(insert_index(end-m+1):insert_index(end-m+2)-1)];
                       end
                    end
                    l = state_input_matrix.*generator_matrix;
                    output_binary = mod(sum(l),2);

                    % e = d.*(int2bit(varargin{2}(1),num_register+1))';
                    % f = d.*(int2bit(varargin{2}(2),num_register+1))';
                    % g = mod([sum(e) sum(f)],2);
                    tre.outputs(i,j) = bi2de(output_binary,'left-msb');     
                end
            end
        % Trellis Structure for Rate 2/3 Feedback Convolutional Encoder
        elseif nargin == 3
            psss;

        end
    end

    function code = convolutional_encoder(u,trellis)
        
        v = log2(trellis.numStates); % number of memory
        k = log2(trellis.numInputSymbols); % number of bits input
        n = log2(trellis.numOutputSymbols); % number of bits output
        initial_states = zeros(1,v);% set initial state to all zero
        u_zero_termination = [u zeros(1,v.*k)]; % zero termination
        %u_zero_termination = u; % zero termination
        %code = zeros(1,length(u_zero_termination)/k*n);% store the final code

        previous_states = bi2de(initial_states,'left-msb')+1;
        
        m = buffer(u_zero_termination, k)';
        input_bits = bi2de(m, 'left-msb')'+1;
        output = zeros(1,length(m));

        for i  = 1:(length(u_zero_termination)/k)
            
            %input_bits = bi2de(u_zero_termination(1+k*(i-1):k*i),'left-msb')+1;
            next_states = trellis.nextStates(previous_states,input_bits(i));
            output(i) = trellis.outputs(previous_states,input_bits(i));
                
            %output_codebits = (int2bit(output,n))';
            %code(1+(i-1)*n:n+(i-1)*n) = output_codebits;

            previous_states = next_states+1;
        
        end
        
        output_codebits = int2bit(output,n);
        code = reshape(output_codebits,1,[]);

    end

    function symb = bits2qpsk(bits)
    %Encode bits as qpsk symbols
        const = [(1+1i),(1-1i),(-1+1i),(-1-1i)]/sqrt(2); % Constellation: QPSK/4-QAM
        m = buffer(bits, 2)';                            % Group bits into bits per symbol
        m_idx = bi2de(m, 'left-msb')'+1;                 % Bits to symbol index
        symb = const(m_idx);                             % Look up symbols using the indices
    end

    function symb = bits2bpsk(bits)
    %Encode bits as bpsk symbols
        bits(bits>0)=-1;    
        bits(bits<=0 & bits>-1)=1;    
        symb = bits;
    end

    function symb = bits2ampm(bits)
    %Encode bits as qpsk symbols
        const = [(1-1i),(-3+3i),(1+3i),(-3-1i),(3-3i),(-1+1i),(3+1i),(-1-3i)]./sqrt(10); 
        m = buffer(bits, 3)';                            % Group bits into bits per symbol
        m_idx = bi2de(m, 'left-msb')'+1;                 % Bits to symbol index
        symb = const(m_idx);                             % Look up symbols using the indices
    end
    
    function y = add_Gaussian_noise(x,sigma)
    %add Gaussian noise to transmitted signal
        if isreal(x)
            noise = sigma*(randn(1,length(x)));
        else 
            noise = sigma*(randn(1,length(x))+1i*randn(1,length(x)));
        end
        y = x + noise;
    end
    
    function c_hat = symbol_detector(y,decodeType) 
        tf = strcmp(decodeType,'hard');
        % bpsk
        if y == real(y)
            % hard detect
            if tf == 1
                y(y>=0) = 0;
                y(y<0) = 1;
            end
            c_hat = y;
        % qpsk
        else
            B = real(y);
            C = imag(y);
            % hard detect
            if tf == 1
                B(B>0)=0;
                B(B<0)=1;
        
                C(C>0)=0;
                C(C<0)=1;
            end
    
            c_hat = zeros(1,2*length(y));
            for i = 1:length(y)
                c_hat(2*i-1) = B(i);
                c_hat(2*i) = C(i);
            end
        end
    end
    
%     function u_hat = hard_input_Viterbi(c_hat,trellis)
% 
%         initial_state = 0;
%         metrics = zeros(1,trellis.numStates);
%         metrics(2:end) = inf;
%         v = log2(trellis.numStates); % number of memory
%         k = log2(trellis.numInputSymbols); % number of bits input
%         n = log2(trellis.numOutputSymbols); % number of bits output
%         L = length(c_hat)/n;
%         path_matrix = zeros(L,trellis.numStates*3);
%         u_hat = zeros(1,L*k);
%         current_state = initial_state;
%         metrics_new = metrics;
%         for i = 2:L        
%             possible_next_states = trellis.nextStates(current_state+1,:);
%             output = (trellis.outputs(current_state+1,:));
%             output_T = output';
%             output2bi = int2bit((output_T(:))',n);
%             received_codebits = (c_hat((1+n*(i-2)):n*(i-1)))';
%             lambda_i_minus_1 = reshape(sum(bitxor(output2bi,received_codebits)), ...
%                                         size(possible_next_states,2), ...
%                                         size(possible_next_states,1))';
%             next_state_value = unique(possible_next_states);
%             %根据每一个可能下一步状态（第i步），求出到它的最小路径，
%             %记录（目的状态，上一步状态和输入bit）到path_matrix
%             for j = 1:length(next_state_value)
%                 %将相同数值的possible_next_states的位置找出来
%                 idx=find(possible_next_states==next_state_value(j));
%                 [row,col] = ind2sub(size(possible_next_states),idx);
%                 %针对同一个目的状态，找到最小的metrics
%                 least_metrics = zeros(1,length(row));
%                 for m = 1:length(row)             
%                     %新的metrics等于它上一个状态的积累metrics+新的汉明距离
%                     least_metrics(m) = metrics(current_state(row(m))+1)+lambda_i_minus_1(row(m),col(m));         
%                 end
%                 %对于每一个j,找到下一步状态的幸存路径，M是指截止下一步状态为止的metrics,I是指它的上一步状态（+1）
%                 [M,I] = min(least_metrics);
%                 metrics_new(next_state_value(j)+1) = M;
%                 %output_bits = output(row(I),col(I));            
%                 %[current_state previous_state inputbit]
%                 path_matrix(L-i+2,3*next_state_value(j)+1:3*next_state_value(j)+3)=[next_state_value(j),current_state(row(I)),col(I)-1];
%             end            
%             current_state = next_state_value;
%             metrics = metrics_new;
%         end
%         
%         %find the survival path from path_matrix
%         target_state = 0;
%         for i = 1:L
%             previou_state = path_matrix(i,3*target_state+1+1);
%             input_bit = path_matrix(i,3*target_state+1+2);
%             u_hat((L-i)*k+1:(L-i)*k+k) = (int2bit(input_bit,k))';
%             target_state = previou_state;
%         end
%         
%         u_hat = u_hat(1:end-v*k);
%                
%     end

    function u_hat = viterbi_decoding(c_hat,trellis,decodeType)

        initial_state = 0;
        metrics = zeros(1,trellis.numStates);
        metrics(2:end) = inf;
        tf = strcmp(decodeType,'hard');
        v = log2(trellis.numStates); % number of memory
        k = log2(trellis.numInputSymbols); % number of bits input
        n = log2(trellis.numOutputSymbols); % number of bits output
        L = length(c_hat)/n;
        % the path_matrix stores the survival path
        % path_matrix = zeros(L,trellis.numStates*3);
        path_matrix = zeros(L,trellis.numStates);
        % u_hat = zeros(1,L*k);
        current_state = initial_state;
        metrics_new = metrics;


        input_bit = zeros(1,L);
        received_codebits = buffer(c_hat, n);
        if tf == 1
            received_codebits2dec = bin2dec(num2str(received_codebits'));

            received_codebits_unique = (unique(received_codebits2dec))';
            received_codebits_unique_2bin = int2bit(received_codebits_unique,n);
        end
        % genarate the lambda_i_minus_1 for each step
        for i = 2:L
            if i < 5
                possible_next_states = trellis.nextStates(current_state+1,:);
                output = (trellis.outputs(current_state+1,:));
                output_T = output';
                output2bi = int2bit((output_T(:))',n);

                % hard decoding hamming distance matrix
                if tf == 1
                    lambda_i_minus_1 = reshape(sum(bitxor(output2bi,received_codebits(:,i-1))), ...
                                                size(possible_next_states,2), ...
                                                size(possible_next_states,1))';
                % soft decoding Euclidean distance matrix
                else
                    minus_expoe = (-1).^output2bi;
                    %minus_expoe = 
                    processed_lambda = sum(minus_expoe.*(received_codebits(:,i-1)));
                    lambda_i_minus_1 = reshape(processed_lambda, ...
                                                size(possible_next_states,2), ...
                                                size(possible_next_states,1))';
                
                end
            
            elseif i == 5
                if tf == 1
                    lambda_matrix = zeros(size(possible_next_states,1),size(possible_next_states,2)*length(received_codebits_unique));
                    for j = 1:length(received_codebits_unique)
    
                        lambda_matrix(:,1+(j-1)*trellis.numInputSymbols:j*trellis.numInputSymbols) = reshape(sum(bitxor(output2bi,received_codebits_unique_2bin(:,j))), ...
                                                                                                                        size(possible_next_states,2), ...
                                                                                                                        size(possible_next_states,1))';
                    end     
                    lambda_i_minus_1 = lambda_matrix(:,1+trellis.numInputSymbols*received_codebits2dec(i-1):trellis.numInputSymbols*(1+received_codebits2dec(i-1)));
                else
                    processed_lambda = sum(minus_expoe.*(received_codebits(:,i-1)));
                    lambda_i_minus_1 = reshape(processed_lambda,size(possible_next_states,2),size(possible_next_states,1))';
                end
                
            else
                if tf == 1  
                    lambda_i_minus_1 = lambda_matrix(:,1+trellis.numInputSymbols*received_codebits2dec(i-1):trellis.numInputSymbols*(1+received_codebits2dec(i-1)));
                else
                    processed_lambda = sum(minus_expoe.*(received_codebits(:,i-1)));
                    lambda_i_minus_1 = reshape(processed_lambda,size(possible_next_states,2),size(possible_next_states,1))';
                end
            end


            next_state_value = unique(possible_next_states);
            
            % According to each state of possible_next_states, find out the
            % survival path, then store (next_states, current_state, and
            % input bit) into path_matrix.

            for j = 1:length(next_state_value)
                % find the location of same value in possible_next_states
                idx=find(possible_next_states==next_state_value(j));
                [row,col] = ind2sub(size(possible_next_states),idx);

                % for the same next state, find the survival path
                % which means find the minimum value of metrics for hard
                % decoding,but find the maximum value of metrics for soft
                % decoding.

                % least_metrics is created to store all possible metrics
                % value for each state for each step
                least_metrics = zeros(1,length(row));      
                for m = 1:length(row) 
                    % new possible metrics for each state equal to previous
                    % metrics of this state + lambda_i_minus_1
                    least_metrics(m) = metrics(current_state(row(m))+1)+lambda_i_minus_1(row(m),col(m));         
                end
                %least_metrics = metrics(current_state(row)+1)+lambda_i_minus_1(row,unique(col)); 
                % for each j, find the survival path of each step
                % M is the next state Metrics, I is the current state(+1)
                if tf == 1
                    [M,I] = min(least_metrics);
                else
                    [M,I] = max(least_metrics);
                end

                metrics_new(next_state_value(j)+1) = M;
                % output_bits = output(row(I),col(I));            
                % [current_state previous_state inputbit]
                % path_matrix(L-i+2,3*next_state_value(j)+1:3*next_state_value(j)+3)=[next_state_value(j),current_state(row(I)),col(I)-1];
                % define : path_matrix = zeros(L,trellis.numStates);
                
                % column value is next state, real part is previous state, imag part is input 
                % path_matrix(L-i+2,I+1) = next_state_value(j) + 1i* (col(I)-1);
                % disp([next_state_value(j)+1 current_state(row(I)) col(I)-1]);

                path_matrix(L-i+2,next_state_value(j)+1) = current_state(row(I)) + 1i* (col(I)-1);
            end            
            current_state = next_state_value;
            metrics = metrics_new;
        end
        
        %find the survival path from path_matrix
        target_state = 0;
        for i = 1:L
            previous_state = real(path_matrix(i,target_state+1));
            input_bit(L+1-i) = imag(path_matrix(i,target_state+1));

            % u_hat((L-i)*k+1:(L-i)*k+k) = (int2bit(input_bit,k))';
            target_state = previous_state;
        end
        
        output_codebits = int2bit(input_bit,k);
        u_hat = reshape(output_codebits,1,[]);

        u_hat = u_hat(1:end-v*k);
               
    end
    
    function u_hat = ampm_detector(y)
        y = y(1:(end-1));
        u_hat = zeros(1,3*length(y)+1);
        a = sqrt(0.1);
        B = real(y);
        C = imag(y);
        for i = 1:length(y)
            if B(i) > C(i)
                if (B(i)+C(i)>(2*a) && C(i)>(-1*a))
                    u_hat(1+3*(i-1):3*i) = [1 1 0];
                elseif (B(i)+C(i)<=(2*a) && C(i)-B(i)>(-4*a) && C(i)+B(i) > -2*a)
                    u_hat(1+3*(i-1):3*i) = [0 0 0];  
                elseif (C(i)<=-1*a && C(i)-B(i)<=(-4*a) && B(i)>=a)
                    u_hat(1+3*(i-1):3*i) = [1 0 0];
                else
                    u_hat(1+3*(i-1):3*i) = [1 1 1];
                end
            else
                if (B(i)+C(i)>(2*a) && B(i)>=-1*a)
                    u_hat(1+3*(i-1):3*i) = [0 1 0]; 
                elseif (C(i)-B(i)>=4*a && B(i)<-1*a && C(i)>=a)
                    u_hat(1+3*(i-1):3*i) = [0 0 1];
                elseif (C(i)<a && C(i)+B(i)<=-2*a)
                    u_hat(1+3*(i-1):3*i) = [0 1 1];
                else
                    u_hat(1+3*(i-1):3*i) = [1 0 1];
                end
            end
        end
    end
    
    function u_hat = viterbi_soft(c_hat,trellis,decodeType)
        initial_state = 0;
        metrics = zeros(1,trellis.numStates);
        metrics(2:end) = inf;
        v = log2(trellis.numStates); % number of memory 3
        k = log2(trellis.numInputSymbols); % number of bits input 2
        n = log2(trellis.numOutputSymbols); % number of bits output 3
        L = length(c_hat);
        % the path_matrix stores the survival path
        path_matrix = zeros(L,trellis.numStates);
        %u_hat = zeros(1,L*k);
        current_state = initial_state;
        metrics_new = metrics;

        input_bit = zeros(1,L);

        tf_ampm = strcmp(decodeType,'ampm');
        tf_bpsk = strcmp(decodeType,'bpsk');
        tf_qpsk = strcmp(decodeType,'qpsk');

        if tf_ampm==1
            const = [(1-1i),(-3+3i),(1+3i),(-3-1i),(3-3i),(-1+1i),(3+1i),(-1-3i)]./sqrt(10);
        end
        if tf_qpsk==1
            const =[(1+1i),(1-1i),(-1+1i),(-1-1i)]/sqrt(2);
        end
        if tf_bpsk==1
            const = [1 -1];
        end
        

        
        for i = 2:L
            if i < 10
                possible_next_states = trellis.nextStates(current_state+1,:);
                output = (trellis.outputs(current_state+1,:));     
                symbol_output = const(output+1);
            end
            received_codebits = c_hat(i-1);
            lambda_i_minus_1 = real(symbol_output-received_codebits).^2 + imag(symbol_output-received_codebits).^2;
            
            next_state_value = unique(possible_next_states);

            for j = 1:length(next_state_value)
                % find the location of same value in possible_next_states
                idx=find(possible_next_states==next_state_value(j));
                [row,col] = ind2sub(size(possible_next_states),idx);

                % for the same next state, find the survival path
                % which means find the minimum value of metrics for hard
                % decoding,but find the maximum value of metrics for soft
                % decoding.

                % least_metrics is created to store all possible metrics
                % value for each state for each step
                least_metrics = zeros(1,length(row));      
                for m = 1:length(row) 
                    % new possible metrics for each state equal to previous
                    % metrics of this state + lambda_i_minus_1
                    least_metrics(m) = metrics(current_state(row(m))+1)+lambda_i_minus_1(row(m),col(m));         
                end

                % for each j, find the survival path of each step
                % M is the next state Metrics, I is the current state(+1)
                [M,I] = min(least_metrics);
                metrics_new(next_state_value(j)+1) = M;
                % output_bits = output(row(I),col(I));            
                % [current_state previous_state inputbit]
                % path_matrix(L-i+2,3*next_state_value(j)+1:3*next_state_value(j)+3)=[next_state_value(j),current_state(row(I)),col(I)-1];
                % define : path_matrix = zeros(L,trellis.numStates);
                
                % column value is next state, real part is previous state, imag part is input 
                % path_matrix(L-i+2,I+1) = next_state_value(j) + 1i* (col(I)-1);
                % disp([next_state_value(j)+1 current_state(row(I)) col(I)-1]);

                path_matrix(L-i+2,next_state_value(j)+1) = current_state(row(I)) + 1i* (col(I)-1);
            end            
            current_state = next_state_value;
            metrics = metrics_new;
        end

        target_state = 0;
        for i = 1:L
            previous_state = real(path_matrix(i,target_state+1));
            input_bit(L+1-i) = imag(path_matrix(i,target_state+1));

            % u_hat((L-i)*k+1:(L-i)*k+k) = (int2bit(input_bit,k))';

            target_state = previous_state;
        end
        
        output_codebits = int2bit(input_bit,k);
        u_hat = reshape(output_codebits,1,[]);

        u_hat = u_hat(1:end-v*k);

    end





% Generate structure with handles to functions
funs.bits2qpsk = @bits2qpsk;
funs.bits2bpsk = @bits2bpsk;
funs.bits2ampm = @bits2ampm;
funs.add_Gaussian_noise = @add_Gaussian_noise;
funs.convolutional_encoder = @convolutional_encoder;
funs.polynomial2trellis = @polynomial2trellis;
funs.symbol_detector = @symbol_detector;
funs.viterbi_decoding = @viterbi_decoding;
funs.ampm_detector = @ampm_detector;
funs.viterbi_soft = @viterbi_soft;

end
function y_rece = viterbi(y)
N = length(y);
nextStates = [0 2;0 2;1 3;1 3];
outputs = [0 3;3 0;1 2;2 1];
dum = 1000;
metric = dum*ones(4,N);             %This is the cumulative matrix the which contains 4 rows corresponding to the four states and N columns corresponding the length of the bit stream ...
% where the real part indicates the minimum path up to the ith element at
% each state and the imaginary part points to the privous node that lead us
% to the current state
state = 1;
for i = 1 : N

    if i == 1
        metric(nextStates(state,1)+1, i) = sum(de2bi(bitxor(y(i),outputs(1,1))))+ 1i*state;
        metric(nextStates(state,2)+1, i) = sum(de2bi(bitxor(y(i),outputs(1,2))))+ 1i*state;
    elseif i == 2
        for j = 1 : 2 : 3
            metric(nextStates(j,1)+1, i) = mod(sum(de2bi(bitxor(y(i),outputs(j,1)))) + real(metric(j,i-1)), dum) + 1i*j;
            metric(nextStates(j,2)+1, i) = mod(sum(de2bi(bitxor(y(i),outputs(j,2)))) + real(metric(j,i-1)), dum) + 1i*j;
        end

    else
        for j = 1: 4
            x= mod(sum(de2bi(bitxor(y(i),outputs(j,1)))) + real(metric(j,i-1)), dum) + 1i*j ; 
            if real(x) < real(metric(nextStates(j,1)+1, i))
                metric(nextStates(j,1)+1, i) = x ; 
            end
            x = mod(sum(de2bi(bitxor(y(i),outputs(j,2)))) + real(metric(j,i-1)), dum) + 1i*j ;
            if real(x) < real(metric(nextStates(j,2)+1, i))
                metric(nextStates(j,2)+1, i) = x ; 
            end
        end
         

    end
end
[~,state] = min(real(metric(:,end)));           % Because our system is not zero-terminated so we take the minimum real value of the Nth column of all states and look at the imaginary number to determine how to go back
y_rece = ones(1,N);
for k = N : -1 : 1
    if state <= 2
        y_rece(k) = 0;
    else
        y_rece(k) = 1;
    end
    state = imag(metric(state,k));
end
end
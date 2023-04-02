% ======================================================================= %
% SSY125 Project
% ======================================================================= %
clc
clear

% ======================================================================= %
% Simulation Options
% ======================================================================= %
N = 1e5;  % simulate N bits each transmission (one block) % In order to have a stable plot we can increase the number of bits
maxNumErrs = 100; % get at least 100 bit errors (more is better)
maxNum = 1e6; % OR stop if maxNum bits have been simulated
EbN0 = -1:8; % power efficiency range

% ======================================================================= %
% Other Options
% ======================================================================= %
% ...

% ======================================================================= %
% Simulation Chain
% ======================================================================= %
BER = zeros(1, length(EbN0)); % pre-allocate a vector for BER results
BER_coded = zeros(1, length(EbN0)); % pre-allocate a vector for BER results

for i = 1:length(EbN0) % use parfor ('help parfor') to parallelize  
  totErr = 0;  % Number of errors observed
  totErr_coded = 0;  % Number of errors observed for the coded
  num = 0; % Number of bits processed
  num_coded = 0; % Number of bits processed for the coded 
  while((totErr < maxNumErrs) && (num < maxNum))
  % ===================================================================== %
  % Begin processing one block of information
  % ===================================================================== %
  % [SRC] generate N information bits 
  bits = randi([0,1],N,1)';
  bits_ref = bits;
  % ... 

  % [ENC] convolutional encoder 
  % Here man can see how trellis is built basically containing two matrises
  %, next-state and output respectively. 

  trellis = struct('nextStates',[0 2;0 2;1 3;1 3],...
  'outputs',[0 3;3 0;1 2;2 1]);
   state = 1;            %Here state = 1 indicates that we start from state 00
   trellis_out = 0;
   Bits_coded = ones(1,2*N); % The encoded bits will be stored in a vector of length 2*N because this generator generates two bits for a signle input bit
   for q = 1:N
       trellis_out = trellis.outputs(state, bits(q) +1);
       state = trellis.nextStates(state, bits(q)+1)+1;
       Bits_coded(2*q -1) = floor(trellis_out/2);
       Bits_coded(2*q) = mod(trellis_out,2);
    
   end



  % ...

  % [MOD] symbol mapper
  bits = Bits_coded;

  constellation = [(1 + 1i), (1 - 1i), (-1 +1i), (-1 -1i)] / sqrt(2);
  BitPar_coded = buffer(bits,log2(length(constellation)))';
  ind = bi2de(BitPar_coded, 'left-msb')'+1;
  symb_coded = constellation(ind);

  BitPar_coded = buffer(bits_ref,log2(length(constellation)))';
  ind = bi2de(BitPar_coded, 'left-msb')'+1;
  symb = constellation(ind);

  % ...

  % [CHA] add Gaussian noise
  Eb = 0.5;                                         %Bit energi for QPSK constellation 
  N0 = Eb./(10.^(EbN0(i)./10));
  N0_coded = 2.*Eb./(10.^(EbN0(i)./10));            % Here we are considering Rc in our calculation because we are using the convolutional encoder of Rc rate
  sigma = sqrt(N0/2);
  sigma_coded = sqrt(N0_coded/2);
  nR = normrnd(0,sigma_coded,[1 length(symb_coded)]);
  nI = normrnd(0,sigma_coded,[1 length(symb_coded)]);
  y_coded = symb_coded+ nR + 1j.*nI;            % Adding Gaussian noise to coded bits 

  nR = normrnd(0,sigma,[1 length(symb)]);
  nI = normrnd(0,sigma,[1 length(symb)]);
  y = symb+ nR + 1j.*nI;                        % Adding Gaussian noise to uncoded bits 
  % ...

  %scatterplot; dont understand why we have scatterplot here or it is just
  %a comment
  plot(y, 'b.')  
  plot(y_coded, 'b.')  

  % [HR] Hard Receiver
  y_rece = hard_receiver(y);
  BitR_coded = hard_receiver(y_coded);
  BitR_de_coded = bi2de(buffer(BitR_coded,2)','left-msb')';
   % Viterbi
  y_rece_coded = viterbi(BitR_de_coded);




   %



  % ...

  % [SR] Soft Receiver
  % ...
  % ===================================================================== %
  % End processing one block of information
  % ===================================================================== %
  BitErrs = sum(y_rece ~= bits_ref); % count the bit errors and evaluate the bit error rate
  totErr = totErr + BitErrs;
  num = num + N; 

  BitErrs_coded = sum(y_rece_coded ~= bits_ref); % count the bit errors and evaluate the bit error rate of the coded stream
  totErr_coded = totErr_coded + BitErrs_coded;
  num_coded = num_coded + N; 


  disp(['+++ ' num2str(totErr) '/' num2str(maxNumErrs) ' errors. '...
      num2str(num) '/' num2str(maxNum) ' bits. Projected error rate = '...
      num2str(totErr/num, '%10.1e') '. +++']);

  disp(['+++ ' num2str(totErr_coded) '/' num2str(maxNumErrs) ' errors. '...
      num2str(num_coded) '/' num2str(maxNum) ' bits. Projected error rate = '...
      num2str(totErr_coded/num, '%10.1e') '. +++']);
  end 
  BER(i) = totErr/num; 
  BER_coded(i) = totErr_coded/num_coded;
end

BERTh = qfunc(sqrt(2*10.^(EbN0./10)));
semilogy(EbN0, BER,'-',EbN0, BERTh, '*',EbN0, BER_coded, 'r')
ylim([1e-4 1]);
xlabel('$E_b / N_0 \ [dB]$','Interpreter','Latex')
ylabel('$BER$','Interpreter','Latex')
legend('Uncoded system (simulation)','Uncoded system (theory)', 'coded system')
% ======================================================================= %
% End
% ======================================================================= %
% ======================================================================= %
% Part2_Coding can Increase Efficiency
% ======================================================================= %
clc
clear

funs = student_sols();

% ======================================================================= %
% Simulation Options
% ======================================================================= %
N = 1e5;  % simulate N bits each transmission (one block)
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

for i = 1:length(EbN0) % use parfor ('help parfor') to parallelize  
  % soft
  totErr = 0;  % Number of errors observed
  num = 0; % Number of bits processed
   

  while((totErr < maxNumErrs) && (num < maxNum))
  % ===================================================================== %
  % Begin processing one block of information
  % ===================================================================== %
  % [SRC] generate N information bits 
    u = randsrc(1,N,[0,1]);

  % [ENC] convolutional encoder    
    trellis4.numInputSymbols = 2^2;
    trellis4.numOutputSymbols = 2^3;
    trellis4.numStates  = 2^3;
    trellis4.nextStates = [0  2  1  3; ...
                           4  6  5  7; ...
                           1  3  0  2; ...
                           5  7  4  6; ...
                           2  0  3  1; ...
                           6  4  7  5; ...
                           3  1  2  0; ...
                           7  5  6  4];
    
    trellis4.outputs =    [0  1  2  3; ...
                           4  5  6  7; ...
                           0  1  2  3; ...
                           4  5  6  7; ...
                           0  1  2  3; ...
                           4  5  6  7; ...
                           0  1  2  3; ...
                           4  5  6  7];


    c4 = funs.convolutional_encoder(u,trellis4);

  % [MOD] symbol mapper

    x_system_3 = funs.bits2ampm(c4);
    

  % [CHA] add Gaussian noise    
    Eb_ampm = 1/3;
    Rc_2 = 2/3;
    sigma_system_3 = sqrt(Eb_ampm/(2*10.^((i-2)/10)));
    y_system_3 = funs.add_Gaussian_noise(x_system_3,sqrt(3/2)*sigma_system_3);
    
    
  % [SR] Soft Receiver
    u_hat_soft_ampm = funs.viterbi_ampm(y_system_3,trellis4);
    

  % ===================================================================== %
  % End processing one block of information
  % ===================================================================== %
  % soft_ampm
  BitErrs = sum(u~=u_hat_soft_ampm); 
  totErr = totErr + BitErrs;
  num = num + N; 
  disp(['+++ ' num2str(totErr) '/' num2str(maxNumErrs) ' errors. '...
      num2str(num) '/' num2str(maxNum) ' bits. Projected error rate = '...
      num2str(totErr/num, '%10.1e') '. +++']);
  

  end 
  BER(i) = totErr/num;

end



figure(3);
% semilogy(EbN0, BER,'-',EbN0, BER_theory, '--',EbN0,BER_uncoded,'*')
semilogy(EbN0,BER)
ylim([1e-4 1]);
xlabel('Eb/N0 [dB]')
ylabel('BER')
%legend('coded system 2','coded system 1','coded system 3','uncoded system 2','uncoded system 1','uncoded system 3')


% ======================================================================= %
% End
% ======================================================================= %
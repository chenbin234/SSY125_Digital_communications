% ======================================================================= %
% Part2_Coding can Increase Efficiency
% ======================================================================= %
clc
clear
tic
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
BER_3_BPSK = zeros(1, length(EbN0));
BER_4_AMPM = zeros(1, length(EbN0));
BER_uncoded_bpsk = zeros(1, length(EbN0));
BER_uncoded_qpsk = zeros(1, length(EbN0));
BER_uncoded_ampm = zeros(1, length(EbN0));

for i = 1:length(EbN0) % use parfor ('help parfor') to parallelize  
  % soft
  totErr = 0;  % Number of errors observed
  num = 0; % Number of bits processed
  totErr_bpsk = 0;
  num_bpsk = 0;
  totErr_ampm = 0;
  num_ampm = 0;
  % hard
  totErr_uncoded_bpsk = 0;
  num_uncoded_bpsk = 0;
  totErr_uncoded_qpsk = 0;
  num_uncoded_qpsk = 0;
  totErr_uncoded_ampm = 0;
  num_uncoded_ampm = 0;    

  while((totErr < maxNumErrs) && (num < maxNum))
  % ===================================================================== %
  % Begin processing one block of information
  % ===================================================================== %
  % [SRC] generate N information bits 
    u = randsrc(1,N,[0,1]);

  % [ENC] convolutional encoder
    trellis3 = funs.polynomial2trellis(5,[19 27]);
    
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


    c3 = funs.convolutional_encoder(u,trellis3);
    c4 = funs.convolutional_encoder(u,trellis4);

  % [MOD] symbol mapper
    x_system_1 = funs.bits2bpsk(c3);
    x_system_2 = funs.bits2qpsk(c3);
    x_system_3 = funs.bits2ampm(c4);
    
    x_uncoded_1 = funs.bits2bpsk(u);
    x_uncoded_2 = funs.bits2qpsk(u);
    x_uncoded_3 = funs.bits2ampm(u);

  % [CHA] add Gaussian noise

    Eb_qpsk = 0.5;
    Rc_1 = 0.5;
    sigma_system_2 = sqrt(Eb_qpsk/(2*10.^((i-2)/10)));
    y_system_2 = funs.add_Gaussian_noise(x_system_2,sqrt(2)*sigma_system_2);
  
    Eb_bpsk = 1;
    Rc_1 = 0.5;
    % sigma_system_1 = sqrt(Eb_bpsk/(2*10.^((i-2)/10)));
    sigma_system_1 = sqrt(Eb_bpsk/(2*10.^((i-2)/10))); % sigma square = N0, since there is no imaginary part;
    y_system_1 = funs.add_Gaussian_noise(x_system_1,sqrt(2)*sigma_system_1);
  
    Eb_ampm = 1/3;
    Rc_2 = 2/3;
    sigma_system_3 = sqrt(Eb_ampm/(2*10.^((i-2)/10)));
    y_system_3 = funs.add_Gaussian_noise(x_system_3,sqrt(3/2)*sigma_system_3);
    
    
    y_uncoded_bpsk = funs.add_Gaussian_noise(x_uncoded_1,sigma_system_1);
    y_uncoded_qpsk = funs.add_Gaussian_noise(x_uncoded_2,sigma_system_2);
    y_uncoded_ampm = funs.add_Gaussian_noise(x_uncoded_3,sigma_system_3);

  % scatterplot: plot(y, 'b.')  
  % plot(y, 'b.') 

  % [HR] Hard Receiver
  % c_hat = funs.symbol_detector(y);
  % u_hat = funs.hard_input_Viterbi(c_hat,trellis);
    u_hat_uncoded_1 = funs.symbol_detector(y_uncoded_bpsk,'hard');
    u_hat_uncoded_2 = funs.symbol_detector(y_uncoded_qpsk,'hard');
    u_hat_uncoded_3 = funs.ampm_detector(y_uncoded_ampm);
    

  % [SR] Soft Receiver
    c_hat_soft_bpsk = funs.symbol_detector(y_system_1,'soft');
    u_hat_soft_bpsk = funs.viterbi_decoding(c_hat_soft_bpsk,trellis3,'soft');

    c_hat_soft_qpsk = funs.symbol_detector(y_system_2,'soft');
    u_hat_soft_qpsk = funs.viterbi_decoding(c_hat_soft_qpsk,trellis3,'soft');

    u_hat_soft_ampm = funs.viterbi_ampm(y_system_3,trellis4);
    

  % ===================================================================== %
  % End processing one block of information
  % ===================================================================== %
  % soft_qpsk
  BitErrs = sum(u~=u_hat_soft_qpsk); 
  totErr = totErr + BitErrs;
  num = num + N; 
  disp(['+++ ' num2str(totErr) '/' num2str(maxNumErrs) ' errors. '...
      num2str(num) '/' num2str(maxNum) ' bits. Projected error rate = '...
      num2str(totErr/num, '%10.1e') '. +++']);
  
  % soft_bpsk
  BitErrs_bpsk = sum(u~=u_hat_soft_bpsk); 
  totErr_bpsk = totErr_bpsk + BitErrs_bpsk;
  num_bpsk = num_bpsk + N; 
  disp(['+++ ' num2str(totErr_bpsk) '/' num2str(maxNumErrs) ' errors. '...
      num2str(num_bpsk) '/' num2str(maxNum) ' bits. Projected error rate = '...
      num2str(totErr_bpsk/num_bpsk, '%10.1e') '. +++']);
  
  % soft_ampm  
  BitErrs_ampm = sum(u~=u_hat_soft_ampm); 
  totErr_ampm = totErr_ampm + BitErrs_ampm;
  num_ampm = num_ampm + N; 
  disp(['+++ ' num2str(totErr_ampm) '/' num2str(maxNumErrs) ' errors. '...
      num2str(num_ampm) '/' num2str(maxNum) ' bits. Projected error rate = '...
      num2str(totErr_ampm/num_ampm, '%10.1e') '. +++']);
  
  % hard_qpsk
  BitErrs_uncoded_qpsk = sum(u~=u_hat_uncoded_2); 
  totErr_uncoded_qpsk = totErr_uncoded_qpsk + BitErrs_uncoded_qpsk;
  num_uncoded_qpsk = num_uncoded_qpsk + N; 
  disp(['+++ ' num2str(totErr_uncoded_qpsk) '/' num2str(maxNumErrs) ' errors. '...
      num2str(num_uncoded_qpsk) '/' num2str(maxNum) ' bits. Projected error rate = '...
      num2str(totErr_uncoded_qpsk/num_uncoded_qpsk, '%10.1e') '. +++']);

  % hard_bpsk
  BitErrs_uncoded_bpsk = sum(u~=u_hat_uncoded_1); 
  totErr_uncoded_bpsk = totErr_uncoded_bpsk + BitErrs_uncoded_bpsk;
  num_uncoded_bpsk = num_uncoded_bpsk + N; 
  disp(['+++ ' num2str(totErr_uncoded_bpsk) '/' num2str(maxNumErrs) ' errors. '...
      num2str(num_uncoded_bpsk) '/' num2str(maxNum) ' bits. Projected error rate = '...
      num2str(totErr_uncoded_bpsk/num_uncoded_bpsk, '%10.1e') '. +++']);
  
  % hard_ampm
  BitErrs_uncoded_ampm = sum(u~=u_hat_uncoded_3); 
  totErr_uncoded_ampm = totErr_uncoded_ampm + BitErrs_uncoded_ampm;
  num_uncoded_ampm = num_uncoded_ampm + N; 
  disp(['+++ ' num2str(totErr_uncoded_ampm) '/' num2str(maxNumErrs) ' errors. '...
      num2str(num_uncoded_ampm) '/' num2str(maxNum) ' bits. Projected error rate = '...
      num2str(totErr_uncoded_ampm/num_uncoded_ampm, '%10.1e') '. +++']);

  end 
  BER(i) = totErr/num;
  BER_3_BPSK(i) = totErr_bpsk/num_bpsk;
  BER_4_AMPM(i) = totErr_ampm/num_ampm;
  BER_uncoded_qpsk(i) = totErr_uncoded_qpsk/num_uncoded_qpsk;
  BER_uncoded_bpsk(i) = totErr_uncoded_bpsk/num_uncoded_bpsk;
  BER_uncoded_ampm(i) = totErr_uncoded_ampm/num_uncoded_ampm;
end



figure(3);
% semilogy(EbN0, BER,'-',EbN0, BER_theory, '--',EbN0,BER_uncoded,'*')
semilogy(EbN0,BER, EbN0,BER_3_BPSK, EbN0,BER_4_AMPM, EbN0,BER_uncoded_qpsk, EbN0,BER_uncoded_bpsk, EbN0,BER_uncoded_ampm)
ylim([1e-4 1]);
xlabel('Eb/N0 [dB]')
ylabel('BER')
legend('coded system 2','coded system 1','coded system 3','uncoded system 2','uncoded system 1','uncoded system 3')

toc
disp(['running time:', num2str(toc)]);
% ======================================================================= %
% End
% ======================================================================= %
% ======================================================================= %
% SSY125 Project
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
BER_uncoded = zeros(1, length(EbN0));

for i = 1:length(EbN0) % use parfor ('help parfor') to parallelize  
  totErr = 0;  % Number of errors observed
  num = 0; % Number of bits processed
  totErr_uncoded = 0;
  num_uncoded = 0;

  while((totErr < maxNumErrs) && (num < maxNum))
  % ===================================================================== %
  % Begin processing one block of information
  % ===================================================================== %
  % [SRC] generate N information bits
    u = randsrc(1,N,[0,1]);

  % [ENC] convolutional encoder
    trellis = funs.polynomial2trellis(3,[5 7]);
    c = funs.convolutional_encoder(u,trellis);

  % [MOD] symbol mapper
  % ...
    x = funs.bits2qpsk(c);
    x_uncoded = funs.bits2qpsk(u);
  % [CHA] add Gaussian noise
    Eb = 0.5;
    Rc = 0.5;
    sigma = sqrt(Eb/(2*10.^((i-2)/10)));
    y = funs.add_Gaussian_noise(x,sqrt(2)*sigma);
    y_uncoded = funs.add_Gaussian_noise(x_uncoded,sigma);
  % scatterplot: plot(y, 'b.')  
  % plot(y, 'b.') 
  % [HR] Hard Receiver
 
    c_hat = funs.symbol_detector(y,'hard');
    u_hat = funs.viterbi_decoding(c_hat,trellis,'hard');
    u_hat_uncoded = funs.symbol_detector(y_uncoded,'hard');

  % [SR] Soft Receiver
  % ...
  % ===================================================================== %
  % End processing one block of information
  % ===================================================================== %
  BitErrs = sum(u~=u_hat); % count the bit errors and evaluate the bit error rate
  totErr = totErr + BitErrs;
  num = num + N; 

  disp(['+++ ' num2str(totErr) '/' num2str(maxNumErrs) ' errors. '...
      num2str(num) '/' num2str(maxNum) ' bits. Projected error rate = '...
      num2str(totErr/num, '%10.1e') '. +++']);
  
  % Bit Error rate for uncoded tramsmission
  BitErrs_uncoded = sum(u~=u_hat_uncoded);
  totErr_uncoded = totErr_uncoded + BitErrs_uncoded;
  num_uncoded = num_uncoded + N;

  end 
  BER(i) = totErr/num;
  BER_uncoded(i) = totErr_uncoded/num_uncoded;
end


BER_theory = qfunc(sqrt(2*10.^(EbN0./10)));



% semilogy(EbN0, BER,'-',EbN0, BER_theory, '--',EbN0,BER_uncoded,'*')
semilogy(EbN0, BER, EbN0, BER_theory,EbN0,BER_uncoded,'*')
ylim([1e-4 1]);
xlabel('Eb/N0 [dB]')
ylabel('BER')
legend('coded system ( simulation )','Uncoded system ( theory )','Uncoded system ( simulation )')
toc
disp(['running time:', num2str(toc)]);

% ======================================================================= %
% End
% ======================================================================= %
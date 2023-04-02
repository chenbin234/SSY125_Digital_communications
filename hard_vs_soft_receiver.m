% ======================================================================= %
% Part2_Hard vs. Soft Receiver
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
BER_soft = zeros(1, length(EbN0));

parfor i = 1:length(EbN0) % use parfor ('help parfor') to parallelize  
  totErr = 0;  % Number of errors observed
  num = 0; % Number of bits processed
  totErr_uncoded = 0;
  num_uncoded = 0;
  totErr_soft = 0;
  num_soft = 0;

  % while((totErr < maxNumErrs) && (num < maxNum))
  while(((totErr < maxNumErrs) && (num < maxNum)) || (totErr_soft < 5))
  % ===================================================================== %
  % Begin processing one block of information
  % ===================================================================== %
  % [SRC] generate N information bits 
    u = randsrc(1,N,[0,1]);

  % [ENC] convolutional encoder
    trellis = funs.polynomial2trellis(5,[23 22]);
    c = funs.convolutional_encoder(u,trellis);

  % [MOD] symbol mapper
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
    c_hat_soft = funs.symbol_detector(y,'soft');
    u_hat_soft = funs.viterbi_decoding(c_hat_soft,trellis,'soft');

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

  disp(['+++ ' num2str(totErr_uncoded) '/' num2str(maxNumErrs) ' errors. '...
      num2str(num_uncoded) '/' num2str(maxNum) ' bits. Projected error rate = '...
      num2str(totErr_uncoded/num_uncoded, '%10.1e') '. +++']);


   % Bit Error rate for uncoded tramsmission
  BitErrs_soft = sum(u~=u_hat_soft);
  totErr_soft = totErr_soft + BitErrs_soft;
  num_soft = num_soft + N;

  disp(['+++ ' num2str(totErr_soft) '/' num2str(maxNumErrs) ' errors. '...
      num2str(num_soft) '/' num2str(maxNum) ' bits. Projected error rate = '...
      num2str(totErr_soft/num_soft, '%10.1e') '. +++']);  


  end 
  BER(i) = totErr/num;
  BER_uncoded(i) = totErr_uncoded/num_uncoded;
  BER_soft(i) = totErr_soft/num_soft;
end



%BER_theory = qfunc(sqrt(2*10.^(EbN0./10)));
% The upper bound of BER_soft
% spect = distspec(trellis,10);
% BER_upper = 0;
% for d=spect.dfree:(spect.dfree+9)
%     Ad = spect.weight(d-spect.dfree+1)*spect.event(d-spect.dfree+1)/N;
%     BER_upper = BER_upper + Ad*qfun(sqrt(2*d*Rc*10.^(EbN0./10)));
% end
trellis = funs.polynomial2trellis(5,[23 22]);
BER_upper = zeros(1, length(EbN0));
Rc = 0.5;
for i = 1:length(EbN0)
    spect = distspec(trellis,20);
    BER_upper_bound = 0;
    for d=spect.dfree:(spect.dfree+19)
        %Ad = spect.weight(d-spect.dfree+1)*spect.event(d-spect.dfree+1)./(100000);
        Ad = spect.weight(d-spect.dfree+1);
%         b = sqrt(2*d*Rc*(10.^(EbN0(i)./10)));
%         a = qfunc(sqrt(2*d*Rc*10.^(EbN0(i)./10)));
        BER_upper_bound = BER_upper_bound + Ad*qfunc(sqrt(2*d*Rc*(10.^(EbN0(i)./10))));
    end
    BER_upper(i) = BER_upper_bound;
end



% semilogy(EbN0, BER,'-',EbN0, BER_theory, '--',EbN0,BER_uncoded,'*')
% semilogy(EbN0, BER, EbN0, BER_theory,EbN0,BER_uncoded,'*')
semilogy(EbN0, BER, EbN0, BER_soft, EbN0, BER_uncoded, EbN0, BER_upper)
ylim([1e-4 1]);
xlabel('Eb/N0 [dB]')
ylabel('BER')
legend('coded system (hard)','coded system (soft)','Uncoded system (simulation)','Upper Bound (soft)')
toc
disp(['running time:', num2str(toc)]);

% ======================================================================= %
% End
% ======================================================================= %
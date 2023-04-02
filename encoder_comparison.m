% ======================================================================= %
% Part2_Encoder_comparison
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
BER_2 = zeros(1, length(EbN0));
BER_3 = zeros(1, length(EbN0));

parfor i = 1:length(EbN0) % use parfor ('help parfor') to parallelize  
  totErr = 0;  % Number of errors observed
  num = 0; % Number of bits processed
  totErr_2 = 0;
  num_2 = 0;
  totErr_3 = 0;
  num_3 = 0;

  % while((totErr < maxNumErrs) && (num < maxNum))
  while(((totErr < maxNumErrs) && (num < maxNum)) || (totErr_2 < 5) || (totErr_3 < 5))
  % ===================================================================== %
  % Begin processing one block of information
  % ===================================================================== %
  % [SRC] generate N information bits 
    u = randsrc(1,N,[0,1]);

  % [ENC] convolutional encoder
    trellis1 = funs.polynomial2trellis(3,[5 7]);
    trellis2 = funs.polynomial2trellis(5,[23 22]);
    trellis3 = funs.polynomial2trellis(5,[19 27]);
    c1 = funs.convolutional_encoder(u,trellis1);
    c2 = funs.convolutional_encoder(u,trellis2);
    c3 = funs.convolutional_encoder(u,trellis3);

  % [MOD] symbol mapper
    x1 = funs.bits2qpsk(c1);
    x2 = funs.bits2qpsk(c2);
    x3 = funs.bits2qpsk(c3);

  % [CHA] add Gaussian noise
    Eb = 0.5;
    Rc = 0.5;
    sigma = sqrt(Eb/(2*10.^((i-2)/10)));
    y1 = funs.add_Gaussian_noise(x1,sqrt(2)*sigma);
    y2 = funs.add_Gaussian_noise(x2,sqrt(2)*sigma);
    y3 = funs.add_Gaussian_noise(x3,sqrt(2)*sigma);

  % scatterplot: plot(y, 'b.')  
  % plot(y, 'b.') 

  % [HR] Hard Receiver
  % c_hat = funs.symbol_detector(y,'hard');
  % u_hat = funs.viterbi_decoding(c_hat,trellis,'hard');
  % u_hat_uncoded = funs.symbol_detector(y_uncoded,'hard');

   
  % [SR] Soft Receiver
%     c_hat_soft_1 = funs.symbol_detector(y1,'soft');
%     c_hat_soft_2 = funs.symbol_detector(y2,'soft');
%     c_hat_soft_3 = funs.symbol_detector(y3,'soft');

%     u_hat_soft_1 = funs.viterbi_decoding(c_hat_soft_1,trellis1,'soft');
%     u_hat_soft_2 = funs.viterbi_decoding(c_hat_soft_2,trellis2,'soft');
%     u_hat_soft_3 = funs.viterbi_decoding(c_hat_soft_3,trellis3,'soft');
    
    u_hat_soft_1 = funs.viterbi_soft(y1,trellis1,'qpsk');
    u_hat_soft_2 = funs.viterbi_soft(y2,trellis2,'qpsk');
    u_hat_soft_3 = funs.viterbi_soft(y3,trellis3,'qpsk');


  % ===================================================================== %
  % End processing one block of information
  % ===================================================================== %
  BitErrs = sum(u~=u_hat_soft_1); % count the bit errors and evaluate the bit error rate
  totErr = totErr + BitErrs;
  num = num + N; 

  disp(['+++ ' num2str(totErr) '/' num2str(maxNumErrs) ' errors. '...
      num2str(num) '/' num2str(maxNum) ' bits. Projected error rate = '...
      num2str(totErr/num, '%10.1e') '. +++']);
  
  % Bit Error rate for uncoded tramsmission
  BitErrs_2 = sum(u~=u_hat_soft_2);
  totErr_2 = totErr_2 + BitErrs_2;
  num_2 = num_2 + N;

  disp(['+++ ' num2str(totErr_2) '/' num2str(maxNumErrs) ' errors. '...
      num2str(num_2) '/' num2str(maxNum) ' bits. Projected error rate = '...
      num2str(totErr_2/num_2, '%10.1e') '. +++']);


   % Bit Error rate for uncoded tramsmission
  BitErrs_3 = sum(u~=u_hat_soft_3);
  totErr_3 = totErr_3 + BitErrs_3;
  num_3 = num_3 + N;

  disp(['+++ ' num2str(totErr_3) '/' num2str(maxNumErrs) ' errors. '...
      num2str(num_3) '/' num2str(maxNum) ' bits. Projected error rate = '...
      num2str(totErr_3/num_3, '%10.1e') '. +++']); 


  end 
  BER(i) = totErr/num;
  BER_2(i) = totErr_2/num_2;
  BER_3(i) = totErr_3/num_3;
end



%BER_theory = qfunc(sqrt(2*10.^(EbN0./10)));
% The upper bound of BER_soft
% spect = distspec(trellis,10);
% BER_upper = 0;
% for d=spect.dfree:(spect.dfree+9)
%     Ad = spect.weight(d-spect.dfree+1)*spect.event(d-spect.dfree+1)/N;
%     BER_upper = BER_upper + Ad*qfun(sqrt(2*d*Rc*10.^(EbN0./10)));
% end
trellis1 = funs.polynomial2trellis(3,[5 7]);
trellis2 = funs.polynomial2trellis(5,[23 22]);
trellis3 = funs.polynomial2trellis(5,[19 27]);

trellis_matrix = [trellis1 trellis2 trellis3];
BER_upper = zeros(3, length(EbN0));
Rc = 0.5;

for j =1:3
    trellis = trellis_matrix(j);
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
        BER_upper(j,i) = BER_upper_bound;
    end
end


% semilogy(EbN0, BER,'-',EbN0, BER_theory, '--',EbN0,BER_uncoded,'*')
% semilogy(EbN0, BER, EbN0, BER_theory,EbN0,BER_uncoded,'*')
semilogy(EbN0, BER, EbN0, BER_2, EbN0, BER_3, EbN0, BER_upper(1,:), EbN0, BER_upper(2,:), EbN0, BER_upper(3,:))
ylim([1e-4 1]);
xlabel('Eb/N0 [dB]')
ylabel('BER')
legend('\bf\epsilon 1','\bf\epsilon 2','\bf\epsilon 3','Upper Bound (\bf\epsilon 1)','Upper Bound (\bf\epsilon 2)','Upper Bound (\bf\epsilon 3)')
toc
disp(['running time:', num2str(toc)]);

% ======================================================================= %
% End
% ======================================================================= %
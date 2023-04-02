funs = student_sols();


trellis1 = funs.polynomial2trellis(3,[5 7]);
trellis2 = funs.polynomial2trellis(5,[23 22]);
trellis3 = funs.polynomial2trellis(5,[19 27]);

trellis_matrix = [trellis1 trellis2 trellis3];
BER_upper = zeros(3, length(EbN0));
Rc = 0.5;

for j =1:3
    trellis = trellis_matrix(j);
    for i = 1:length(EbN0)
        spect = distspec(trellis,10);
        BER_upper_bound = 0;
        for d=spect.dfree:(spect.dfree+9)
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
semilogy(EbN0, BER, EbN0, BER_2, EbN0, BER_3, EbN0, BER_upper(1,:), '--',EbN0, BER_upper(2,:),'--', EbN0, BER_upper(3,:),'--')
ylim([1e-4 1]);
xlabel('Eb/N0 [dB]')
ylabel('BER')
legend('\bf\epsilon 1','\bf\epsilon 2','\bf\epsilon 3','Upper Bound (\bf\epsilon 1)','Upper Bound (\bf\epsilon 2)','Upper Bound (\bf\epsilon 3)')



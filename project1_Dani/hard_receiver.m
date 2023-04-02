function BitR = hard_receiver(y)
BitR = zeros(1,2*length(y));
BitR(1:2:end) = (real(y) < 0)*1;
BitR(2:2:end) = (imag(y) < 0)*1;
end
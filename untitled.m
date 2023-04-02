funs = student_sols();

trellis1 = funs.polynomial2trellis(3,[5 7]);
trellis2 = funs.polynomial2trellis(5,[23 22]);
trellis3 = funs.polynomial2trellis(5,[19 27]);


spect1 = distspec(trellis1,5);
d_min1 = spect1.dfree

spect2 = distspec(trellis2,5);
d_min2 = spect2.dfree

spect3 = distspec(trellis3,5);
d_min3 = spect3.dfree


gain1 = 10*log10(0.5*d_min1)

gain2 = 10*log10(0.5*d_min2)

gain3 = 10*log10(0.5*d_min3)
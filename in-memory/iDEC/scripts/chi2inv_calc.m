% % This MATLAB script helps calculate inverse CDF of Chi Sqaured distribution
fp = fopen('chi2inv.txt', 'w');

for delta = [0.9 0.95 0.99]
    fprintf(fp, 'm,delta,x\n');

    for m = 4:2:16
        r = chi2inv(delta, m);
        fprintf(fp, '%d,%.2f,%.6f\n', m, delta, r);
    end

end

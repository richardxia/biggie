% PLOTLETORRESULTS Plots the results of the Letor experiments

data = load('accuracies.txt');

partitions = data(:,1);
false_pos = data(:,2);
false_neg_high = data(:,4);
false_neg_low = data(:,5);

linewidth = 2
h_pos = loglog(partitions, false_pos, 'ko-');
hold on;
h_pos = plot(partitions, false_pos, 'ko-' );
h_neg_high = plot(partitions, false_neg_high, 'bs-');
h_neg_low = plot(partitions, false_neg_low, 'r^-');
set(h_pos, 'linewidth', linewidth);
set(h_neg_high, 'linewidth', linewidth);
set(h_neg_low, 'linewidth', linewidth);

set(gca, 'fontname', 'times', 'fontsize', 18);
h = xlabel('Weirdness threshold');
set(h, 'fontname', 'times', 'fontsize', 20);
h = ylabel('Number of mismatch');
set(h, 'fontname', 'times', 'fontsize', 20);
legend([h_pos, h_neg_high, h_neg_low], 'False-positives', 'False-negatives in high-complexity region', 'False-negatives in low-complexity region', 'location', 'northeast');
print -depsc2 AccuracyCompare.eps;

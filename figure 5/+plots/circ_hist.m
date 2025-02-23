function [mu, std, pval, z, h] = circ_hist(phases, color, fontSize)

h = polarhistogram(phases, 20, "FaceColor", color, "EdgeColor", color, "Normalization", "probability");

mu = circ_mean(phases,[]);
std = circ_std(phases, []);

hold on;
polarplot([mu mu], rlim, "LineWidth", 2, "Color", "b")

[pval, z] = circ_rtest(phases);

ax = gca;
ax.ThetaTick = [0 90 180 270];
ax.ThetaTickLabel = ["0^\circ", "90^\circ", "\pm180^\circ", "-90^\circ"];
ax.FontSize = fontSize;

lim = ax.RLim;
ax.RTick = lim(2);
ax.RTickLabel = string(lim(2));

end


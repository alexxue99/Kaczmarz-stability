fig=figure;
hold on;
h1 = plot(nan, nan, 'r', 'LineWidth', 4);
h2 = plot(nan, nan, 'b', 'LineWidth', 4);
h3 = plot(nan, nan, 'g', 'LineWidth', 4);
h4 = plot(nan, nan, 'black', 'LineWidth', 4);
h5 = plot(nan, nan, 'Color', [0 100/256 50/256]);
legend([h1 h2 h3 h4 h5], {'RK', sprintf('RK + IR'), 'RAK', 'RAK + IR', sprintf('A\\b')}, ...
'Location', 'northoutside', 'Orientation', 'vertical');
axis off;
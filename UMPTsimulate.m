%estimate the true distribution
n = 10;
for i = 1:10000
    x = frnd(2,60,[n 1]);
    sumx(i) = sum(x);
end

figure;histogram(2*sumx,'Normalization','cdf');

x = 0:0.1:(max(2*sumx));
y = chi2cdf(x,2*n);
hold on;plot(x,y)

figure;clear g
g = gramm('x',2*sumx);
g.stat_density("function","cdf");
g.draw();

g.update('x',x,'y',y);
g.geom_line();
g.draw();
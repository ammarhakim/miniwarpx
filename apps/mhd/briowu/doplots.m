X  = load('bw_wv/frame.x');
fd = load('bw_wv/frame.q5');

pr = pressure(5./3.,fd);


__gnuplot_set__ terminal png
__gnuplot_set__ output "bw.png"

subplot(4,2,1); plot(X,fd(:,1),';rho;');
subplot(4,2,2); plot(X,fd(:,2)./fd(:,1),';u;');
subplot(4,2,3); plot(X,fd(:,3)./fd(:,1),';v;');
subplot(4,2,4); plot(X,fd(:,4)./fd(:,1),';w;');
subplot(4,2,5); plot(X,pr,';pr;');
subplot(4,2,6); plot(X,fd(:,6),';bx;');
subplot(4,2,7); plot(X,fd(:,7),';by;');
subplot(4,2,8); plot(X,fd(:,8),';bz;');

oneplot
__gnuplot_set__ terminal x11
__gnuplot_set__ output

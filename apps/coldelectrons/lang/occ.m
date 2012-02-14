function [ux,ne] = occ(ns)

ux = zeros(1,ns+1);
ne = zeros(1,ns+1);
for i = 0:ns
  fn = sprintf('lang/frame.q%d', i);
  fd = load(fn);
  ux(i+1) = fd(20,2)/fd(20,1);
  ne(i+1) = fd(20,1);
end

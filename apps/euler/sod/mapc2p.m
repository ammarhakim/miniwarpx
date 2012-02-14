function Xp = mapc2p(Xc)

nx = max(size(Xc));
Xp = zeros(1,nx);

y0 = 0.5;
r = 5;
h = Xc(end);

B = log( (1+(exp(r)-1)*y0/h) / (1-(1-exp(-r))*y0/h) )/(2*r);

for ix = 1:nx
  eta = Xc(ix);
  Xp(ix) = y0*(1 + sinh(r*(eta-B))/sinh(r*B));
end


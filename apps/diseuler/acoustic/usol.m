function ue = usol(X,cs,wc,t,nt)

nn = max(size(X));
ue = zeros(1,nn);

for n = 0:nt
  kn = 2*pi*(2*n+1); % for a "square" pulse
  wn = sqrt(kn^2*cs^2+wc^2);
  ue = ue + 1/(2*n+1)*sin(kn*X + wn*t);
end

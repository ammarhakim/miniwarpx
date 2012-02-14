function pr = pressure(gas_gamma, fd)

u2 = (fd(:,2).^2 + fd(:,3).^2 + fd(:,3).^2)./fd(:,1);
b2 = fd(:,6).^2 + fd(:,7).^2 + fd(:,8).^2;

pr = (gas_gamma-1)*(fd(:,5) - 0.5*u2 - 0.5*b2);

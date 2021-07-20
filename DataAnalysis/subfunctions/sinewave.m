function y=sinewave(beta,x)

y = beta(1)+(beta(2)-beta(1)).*cos(x+beta(4));
% y=beta(1)+(beta(2)-beta(1)).*cos(2*pi*beta(4)*x+beta(3));


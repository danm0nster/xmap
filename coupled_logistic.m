function [ x, y] = coupled_logistic( x1, y1, rx, ry, bxy, byx, N )
%COUPLED_LOGISTIC(x1, y1, rx, ry, bxy, byx, N) Computes N iteration of the
% coupled logistic map, with initial values x1, y1; control parameters
% (growth rates) rx, ry and coupling between the variables bxy and byx.
x = zeros(1,N);
y = zeros(1,N);
x(1) = x1;
y(1) = y1;
for t = 2:N
    x(t) = x(t-1)*(rx-rx*x(t-1)-bxy*y(t-1));
    y(t) = y(t-1)*(ry-ry*y(t-1)-byx*x(t-1));
end

end


function x = twomoons_model(bigtheta,~,numsim)

if nargin == 1
    numsim=1;
end

if nargin < 2
    numsim=1;
end

t1 = bigtheta(1);
t2 = bigtheta(2);

a = unifrnd(-pi/2,pi/2);
r = 0.1 + 0.01*randn;
p = [r*cos(a)+0.25, r*sin(a)];

x = p + [-abs(t1+t2)/sqrt(2), (-t1+t2)/sqrt(2)];

x = x(:);

end


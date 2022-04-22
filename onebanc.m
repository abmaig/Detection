clear all; close all; clc
syms A1 A2 B1 B2 C1 C2 D1 D2 y h k E I c d J z a b s m

left_bc=0;                      % h = 0 (clamped end)
right_bc=1;                     % h = 1 (free end)
dam_loc=0.5;                    % h = 0.5 (damage location)

d = 10;                         % d = Beam depth (in mm)
p = 20;                          % p = damage percentage = [1-(Id/Iu)]*100
a = 10*(1-(1-(p/100))^(1/3));   % a = crack depth (in mm)
z = a/d;                        % z = Crack depth ratio
E = 200*10^3;                   % E = Elastic modulus (in MPa)
b = 50;                         % b = Width of beam (in mm)

I = (b*(d-a)^3)/12;             % I = Moment of inertia
J = 1.8624*(z^2)-3.95*(z^3)+16.37*(z^4)-37.226*(z^5)+76.81*(z^6)-126.9*(z^7)+172*(z^8)-43.97*(z^9)+66.56*(z^10);
c = 5.346*d*J/(E*I);
k = (1/c);                      % k = Bending constant of the spring

w = [cos(y*h) sin(y*h) cosh(y*h) sinh(y*h)];
% wn = [sin(y*h) cosh(y*h) sinh(y*h)];
XX = [A1;B1;C1;D1;A2;B2;C2;D2];
w1 = w*XX(1:4);
w2 = w*XX(5:8);
% w1 = A1*cos(y*h) + B1*sin(y*h) + C1*cosh(y*h) + D1*sinh(y*h);
% w2 = A2*cos(y*h) + B2*sin(y*h) + C2*cosh(y*h) + D2*sinh(y*h);
% k=1.310616*10^8/eps;
% EI=21.8453*10^10;
% h = x/L
% x=crack distance from clamped end
% L=beam length (500 mm)

d1 = diff(w1,h)/y;                % d = 1st order derivative
d2 = diff(w2,h)/y;
dd1 = diff(w1,h,2)/y^2;           % dd =  2nd order derivative
dd2 = diff(w2,h,2)/y^2;
ddd1 = diff(w1,h,3)/y^3;          % ddd = 3rd order derivative
ddd2 = diff(w2,h,3)/y^3;

% boundary conditions
eqn1 = subs(w1,h,0) == 0;
eqn2 = subs(d1,h,0) == 0;
eqn3 = subs(dd2,h,1) == 0;
eqn4 = subs(ddd2,h,1) == 0;

% compatibility conditions at the location of crack
eqn5 = subs(w1,h,0.5)-subs(w2,h,0.5)==0;
eqn6 = (k*(subs(d1,h,0.5)-subs(d2,h,0.5))+E*I*subs(dd1,h,0.5))==0;
eqn7 = (subs(dd1,h,0.5)-subs(dd2,h,0.5))==0;
eqn8 = (subs(ddd1,h,0.5)-subs(ddd2,h,0.5))==0;

eqns=[eqn1 eqn2 eqn3 eqn4 eqn5 eqn6 eqn7 eqn8];
[D,C] = equationsToMatrix(eqns, XX);
% vars = symvar(eqns);

% X = linsolve(D,C);

%    for yy=0:0.1:10
%       z=subs(det(D),y,yy);
%    end

%   plot(y,z)


f = det(D);
% (subs(f,y,1.4751803509612516141940893578884))

%  y=1;
%
for n = 1:12
    %     s=vpasolve(f==0,y,'random',true)
    s(n) = vpasolve(f==0,y,n);
    m(n) = subs(f,y,s(n));
end

% soly = solve(f==0, y)
% [soly, param, cond] = solve(f==0, y, 'ReturnConditions', true)
% solve(taylor(f,y,'order',100),y,3)
% solve(f==0, y, 'IgnoreAnalyticConstraints', true)
% assume(cond)
% solk = solve(1<soly, soly<20, param)

s_value = s(2);                         % y value (y = [(omega^2)*ro*A/(E*I))^1/4]*L
% L = 500; % in mm
% rho = 800;
% omega = ((s_value/L)^2)*((E*I/rho*A)^(1/2));
subs(f,y,s_value)
plot(s)

D_value = subs(D,y,s_value);
det(D_value)
% rank(D_value)

% elimination approach
D_value(1,:)=[];                         % remove 1st row of D_value
D_colval=D_value(:,1);                   % 1st column of D_value (after removing 1st row)
D_value(:,1)=[];                         % remove 1st column  of D_value
C(1,:)=[];                               % remove 1st row of 'C' column vector
C_n(:,1)=C(:,1)-(D_colval*1);            % subtract 1st column of D_value (after removing 1st row)in 'C' column vector

XX_value=inv(D_value)*C_n;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         XX_value=D_value\C_n;
XX_value_n = double([1;XX_value]);               % taken A1 as unity


w_value = subs(w,y,s_value);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               w_value = subs(w,y,s_value);

w1_value=w_value*XX_value_n(1:4)/2;  % w1_value=subs(w1_value,h,dam_loc);
w2_value=w_value*XX_value_n(5:8)/2;  % w2_value=subs(w2_value,h,dam_loc);
dx = 0.005;
span_L=left_bc:dx:dam_loc;
w1_value_m=(double(subs(w1_value,h,span_L)))';
% w1_value=todecimal(w1_value);
% noise_w1_m = awgn(w1_value_m,SNR);
% figure(111);plot(span_L,w1_value_m,span_L,noise_w1_m);
% legend('original','noise')
% hold on
span_R=dam_loc:dx:right_bc;
w2_value_m=(double(subs(w2_value,h,span_R)))';
% w2_value=todecimal(w2_value);
% figure(111);plot(span_R,w2_value_m);
span = [span_L span_R];
span(:,101) = [];
mode_shape = [w1_value_m;w2_value_m];
mode_shape(101,:) = [];
SNR = 80;
noise_mode_shape = awgn(mode_shape,SNR);
subplot(3,1,1)
plot(span,noise_mode_shape);
legend('original','noise')

% xn = mode_shape + 0.05*randn(1,1);
% plot(mode_shape);
% hold on
% plot(xn,'r')
% xn = mode_shape + 0.5*randn(1,10);
% plot(xn,'r');

% curvature mode shape
curv_mode_shape = diff(noise_mode_shape,2);
curv_mode_shape = [curv_mode_shape;0;0];
subplot(3,1,2)
plot(span,curv_mode_shape);
legend('curvature mode shape')
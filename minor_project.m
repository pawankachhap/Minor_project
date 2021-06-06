%% MATLAB program to determine power absorption by biological bodies
clc
clear all
close all;
a = 2e-2;
an = a*(0.75/pi)^(1/3);

eo = 8.854e-12;
muo = pi*4e-7;
f = 6e9;
er = 80;
sig = 0.84;
w = 2*pi*f;
ko = w*sqrt(muo*eo);
N = 41;

%% The coordinates of the cells are defined here

xd = [0.5*a*ones(1,41)];

yd = [(3*a:-a:-4*a) (4*a:-a:-4*a) (4*a:-a:-4*a) (4*a:-a:-4*a) (1*a:-a:-4*a)];

zd = [-2*a*ones(1,8) -1*a*ones(1,9) 0*a*ones(1,9) 1*a*ones(1,9) 2*a*ones(1,6)];

r = [xd(:) yd(:) zd(:)];

%% Other parameters

tow = [(sig+1i*w*eo*(er-1))*ones(41,1)];

A = -1i*w*muo*ko*a^3/4/pi ;
C = -2*1i*w*muo/3/ko^2*(exp(-1i*ko*an)*(1+1i*ko*an)-1)-1/3/1i/w/eo ;

Ei = [0*exp(-1i*ko*xd(:));86.83*exp(-1i*ko*xd(:));0*exp(-1i*ko*xd(:))];

%%Gxx
for m =1:N
	for n=1:N
		if m == n
			Gxx(m,n)= C*tow(n)-1;
		else
			Rmn= norm(r(m,:)-r(n,:));
			almn = ko*Rmn;
			cxpmn = (r(m,1)-r(n,1))/Rmn;
			cxqmn = (r(m,1)-r(n,1))/Rmn;
			B = (almn^2-1-1i*almn)+cxpmn*cxqmn*(3-almn^2+3*1i*almn);
			Gxx(m,n) = A*tow(n)*exp(-1i*almn)/almn^3*B;
		end
	end
end

%%Gyy
for m =1:N
	for n=1:N
          	if m == n
				  Gyy(m,n)= C*tow(n)-1;
			  else
				Rmn= norm(r(m,:)-r(n,:));
			    almn = ko*Rmn;
				cxpmn = (r(m,2)-r(n,2))/Rmn;
				cxqmn = (r(m,2)-r(n,2))/Rmn;
				B = (almn^2-1-1i*almn)+cxpmn*cxqmn*(3-almn^2+3*1i*almn);
				Gyy(m,n) = A*tow(n)*exp(-1i*almn)/almn^3*B;
			end
	end
end

%%Gzz
for m =1:N
	for n=1:N
        if m == n
			Gzz(m,n)= C*tow(n)-1;
		else
			Rmn= norm(r(m,:)-r(n,:));
			almn = ko*Rmn;
			cxpmn = (r(m,3)-r(n,3))/Rmn;
			cxqmn = (r(m,3)-r(n,3))/Rmn;
			B = (almn^2-1-1i*almn)+cxpmn*cxqmn*(3-almn^2+3*1i*almn);
			Gzz(m,n) = A*tow(n)*exp(-1i*almn)/almn^3*B;
		end
	end
end

%%Gxy
for m=1:N
	for n=1:N
		if m~=n
			Rmn= norm(r(m,:)-r(n,:));
			almn = ko*Rmn;
			cxpmn = (r(m,1)-r(n,1))/Rmn;
			cxqmn = (r(m,2)-r(n,2))/Rmn;

			B = cxpmn*cxqmn*(3-almn^2+3*1i*almn);
			Gxy(m,n) = A*tow(n)*exp(-1i*almn)/almn^3*B;
		end
	end
end

%%Gxz
for m=1:N
	for n=1:N
		if m~=n
			Rmn= norm(r(m,:)-r(n,:));
			almn = ko*Rmn;
			cxpmn = (r(m,1)-r(n,1))/Rmn;
			cxqmn = (r(m,3)-r(n,3))/Rmn;

			B = cxpmn*cxqmn*(3-almn^2+3*1i*almn);
			Gxz(m,n) = A*tow(n)*exp(-1i*almn)/almn^3*B;
		end
	end
end

%Gyx
for m=1:N
	for n=1:N
		if m~=n
			Rmn= norm(r(m,:)-r(n,:));
			almn = ko*Rmn;
			cxpmn = (r(m,2)-r(n,2))/Rmn;
			cxqmn = (r(m,1)-r(n,1))/Rmn;

			B = cxpmn*cxqmn*(3-almn^2+3*1i*almn);
			Gyx(m,n) = A*tow(n)*exp(-1i*almn)/almn^3*B;
		end
	end
end

%%Gyz
for m=1:N
	for n=1:N
		if m~=n
			Rmn= norm(r(m,:)-r(n,:));
			almn = ko*Rmn;
			cxpmn = (r(m,2)-r(n,2))/Rmn;
			cxqmn = (r(m,3)-r(n,3))/Rmn;

			B = cxpmn*cxqmn*(3-almn^2+3*1i*almn);
			Gyz(m,n) = A*tow(n)*exp(-1i*almn)/almn^3*B;
		end
	end
end

%Gzx
for m=1:N
	for n=1:N
		if m~=n
			Rmn= norm(r(m,:)-r(n,:));
			almn = ko*Rmn;
			cxpmn = (r(m,3)-r(n,3))/Rmn;
			cxqmn = (r(m,1)-r(n,1))/Rmn;

			B = cxpmn*cxqmn*(3-almn^2+3*1i*almn);
			Gzx(m,n) = A*tow(n)*exp(-1i*almn)/almn^3*B;
		end
	end
end

%%Gzy
for m=1:N
	for n=1:N
		if m~=n
			Rmn= norm(r(m,:)-r(n,:));
			almn = ko*Rmn;
			cxpmn = (r(m,3)-r(n,3))/Rmn;
			cxqmn = (r(m,2)-r(n,2))/Rmn;

			B = cxpmn*cxqmn*(3-almn^2+3*1i*almn);
			Gzy(m,n) = A*tow(n)*exp(-1i*almn)/almn^3*B;
		end
	end
end

G = [Gxx Gxy Gxz; Gyx Gyy Gyz; Gzx Gzy Gzz];
E = -G\Ei;

%% PLOTS

num = [1:1:41];
%plot(num,abs(Ei));
Ey = E(42:82);
grid on;
plot(num,abs(Ey), 'linewidth',1.5);
xlim([1 41]);
xlabel('Cell number');
ylabel('E(V/m)');

hold on;
Ez = E(83:123);
plot(num,abs(Ez), '--', 'linewidth', 1.5);

legend('y-component','z-component');
title('Absorbed Electric field (Frequency = 6 GHz)');

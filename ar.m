function [ll, lls, j]= ar(parameters,data,backCast)
% data=readtable("C:\Users\thoma\OneDrive\Documents\Bachelorprojekt\per20minut.xlsx")

% logreturn=diff(log(data.Close))
T=size(data,1);
j=zeros(T,1);
a = parameters(1) ;
b = parameters(2) ;
 sigma2 = parameters(3) ;
% d=parameters(3);
% e=parameters(3);
% fejlled=data;
% backCast=0
% parameter=[parameter_0 parameter_1]
% j(1)=data;
j=data;
for i=2:T
    fejlled(i)=j(i)-a*j(i-1)-b;
end
v=8100;
lls=0.5*(log(pi*2)+log(sigma2)+(fejlled.^2)./sigma2);

% lls=0.5*(log(pi*2)+log(sigma2)+(fejlled.^2)./sigma2);
% lls=0.5*(log(sigma2)+(v+1)*log(1+fejlled.^2./(sigma2.*(v-1))));

ll=sum(lls);
end

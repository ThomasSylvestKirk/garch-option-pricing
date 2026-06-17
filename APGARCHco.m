function[ll,lls,h,fejlled]=APGARCH(parameters,logreturn,backCast)
delta=parameters(1);
alpha=parameters(2);
gamma=parameters(3);
omega=parameters(4);
beta=parameters(5);
lambda1=parameters(6);
T=size(logreturn,1);
h=zeros(T,1);
x=zeros(T,1);
h(1)=backCast;
% fejlled=logreturn-e+h/2;
fejlled=logreturn-0.003/251+h/2-lambda1*sqrt(h);

v=length(logreturn)-length(parameters);
% for i =2:T
%     h(i)=(omega + alpha*(abs(fejlled(i-1)-lambda1*sqrt(h(i-1)))- gamma*(fejlled(i-1)-lambda1*sqrt(h(i-1))))^delta+beta*h(i-1)^(delta/2))^(2/delta);
% end
for i =2:T
    h(i)=(omega + alpha*(abs(fejlled(i-1))- gamma*(fejlled(i-1)))^delta+beta*h(i-1)^(delta/2))^(2/delta);
end
lls= 0.5*(log(pi*2)+log(h)+fejlled.^2./h);
% if v<50
%     v = 2 + exp(v);
% else
%     v = 50 + log(v);
% end

% lls=0.5*(log(h)+(v+1)*log(1+fejlled.^2./(h.*(v-1))));

ll= sum(lls);
end



function H= hessian_beggesider(funktion,parametre,varargin)
if size(parametre,2)>size(parametre,1)
     parametre=parametre';
end

f=feval(funktion,parametre,varargin{:});
s=abs(parametre)*0.00001;
parametreh=s+parametre;
s=parametreh-parametre;
s=diag(s);

antal=size(parametre);

afledtep=zeros(antal);
afledten=zeros(antal);
for i=1:antal
afledtep(i)=feval(funktion, parametre+s(:,i),varargin{:});
afledten(i)=feval(funktion, parametre-s(:,i),varargin{:});
end
anafledtep=zeros(antal);
anafledten=zeros(antal);

for j=1:antal
    for i=j:antal
    anafledtep(i,j)=funktion(parametre+s(:,j)+s(:,i),varargin{:});
    anafledtep(j,i)=anafledtep(i,j);
    anafledten(i,j)=funktion(parametre-s(:,j)-s(:,i),varargin{:});
    anafledten(j,i)=anafledten(i,j);
    end
end
H=zeros(antal);
diagonal=diag(s);
GH=diagonal*diagonal';
for i=1:antal
    for j=i:antal
        H(i,j)=(anafledtep(i,j)+anafledten(i,j)-afledtep(i)-afledten(i)-afledtep(j)-afledten(j)+f+f)/GH(i,j)/2;
        H(j,i)=H(i,j);
    end
end
end



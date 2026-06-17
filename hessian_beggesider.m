function H= hessian_beggesider(funktion,parametre,varargin)
if size(parametre,2)>size(parametre,1)
     parametre=parametre';
end
% score=feval(funktion,parametre,logreturn);
% skift=0.0005*estimates;
f=feval(funktion,parametre,varargin{:});
h=abs(parametre)*0.00001;
parametreh=h+parametre;
h=parametreh-parametre;
h=diag(h);

antal=size(parametre);

afledtep=zeros(antal)
afledten=zeros(antal)
for i=1:antal
afledtep(i)=feval(funktion, parametre+h(:,i),varargin{:});
afledten(i)=feval(funktion, parametre-h(:,i),varargin{:});
end
anafledtep=zeros(antal);
anafledten=zeros(antal);

for j=1:antal
    for i=j:antal
    anafledtep(i,j)=funktion(parametre+h(:,j)+h(:,i),varargin{:});
    anafledtep(j,i)=anafledtep(i,j);
    anafledten(i,j)=funktion(parametre-h(:,j)-h(:,i),varargin{:});
    anafledten(j,i)=anafledten(i,j);
    end
end
H=zeros(antal);
diagonal=diag(h);
GH=diagonal*diagonal';
for i=1:antal
    for j=i:antal
        H(i,j)=(anafledtep(i,j)+anafledten(i,j)-afledtep(i)-afledten(i)-afledtep(j)-afledten(j)+f+f)/GH(i,j)/2;
        H(j,i)=H(i,j);
    end
end
end

function [c,ceq]=nlc(parameters,logreturn,fejlled,backCast )
c=parameters(2)*nanmean(((abs(fejlled)).^parameters(1)))./2*((1+parameters(3)).^parameters(1)+(1-parameters(3).^parameters(1)))+parameters(5)-1;
ceq=[];
end
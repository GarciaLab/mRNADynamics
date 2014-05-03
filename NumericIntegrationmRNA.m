
function f = NumericIntegrationmRNA(t,y,halflife,ElongationTime,alpha,Time)

lambda=log(2)/halflife;

if t-ElongationTime/2<0
alphaI=0;
else
alphaI = interp1(Time,alpha,t-ElongationTime/2);
end

f=alphaI-y*lambda;

end


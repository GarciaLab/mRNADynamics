function y=MultipleLinearStepGen(x,X,m,c)

y=zeros(length(x),1);

for i=2:length(x)-1;
    
    index = find(x(i)./X>1,1,'last');
    
    y(i)=m(index)*(x(i))+c(index);
    
    II(i)=index;
end

    save_to_base(1)
    
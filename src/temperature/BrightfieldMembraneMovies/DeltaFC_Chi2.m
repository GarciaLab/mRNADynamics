function Chi2 = DeltaFC_Chi2(x, y,MatSize,LengthY)
%%

TShifts = [0, floor(x(1:MatSize(1)-1))]-min([0, floor(x(1:MatSize(1)-1))])+1;
DeltaShifts =[0, x(MatSize(1):end)];
ymat = NaN(MatSize);
for i = 1:MatSize(1)
    ymat(i,TShifts(i):TShifts(i)+LengthY-1) = y(i,:)+DeltaShifts(i);
end
meanvector =mean(ymat,1,'omitnan');
Chi2 = sum((ymat-meanvector).^2, 'all','omitnan');
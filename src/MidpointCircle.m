% Draw a circle in a matrix using the integer midpoint circle algorithm
% Does not miss or repeat pixels
% Created by : Peter Bone
% Created : 19th March 2007
function i = MidpointCircle(i, radius, xc, yc, value)

%HG modified this function so that it wouldn't freak out if the circle
%touches the edges of the image. It also fills the circle

[Rows,Cols]=size(i);

xc = int16(xc);
yc = int16(yc);

x = int16(0);
y = int16(radius);
d = int16(1 - radius);

if (yc+y<Cols)
    i(xc, yc+y) = value;
else
    i(:,Cols)=value;
end
if (yc-y>0)
    i(xc, yc-y) = value;
else
    i(:,1)=value;
end
if (xc+y<Rows)
    i(xc+y, yc) = value;
else
    i(Rows,:)=value;
end
if (xc-y>0)
    i(xc-y, yc) = value;
else
    i(1,:)=value;
end


while ( x < y - 1 )
    x = x + 1;
    if ( d < 0 ) 
        d = d + x + x + 1;
    else 
        y = y - 1;
        a = x - y + 1;
        d = d + a + a;
    end
    if (x+xc<Rows)&(y+yc<Cols)
        i( x+xc,  y+yc) = value;
    end
    if (y+xc<Rows)&(x+yc<Cols)
        i( y+xc,  x+yc) = value;
    end
    if (y+xc<Rows)&(-x+yc>0)
        i( y+xc, -x+yc) = value;
    end
    if (x+xc<Rows)&(-y+yc>0)
        i( x+xc, -y+yc) = value;
    end
    if (-x+xc>0)&(-y+yc>0)
        i(-x+xc, -y+yc) = value;
    end
    if (-y+xc>0)&(-x+yc>0)
        i(-y+xc, -x+yc) = value;
    end
    if (-y+xc>0)&(x+yc<Cols)
        i(-y+xc,  x+yc) = value;
    end
    if (-x+xc>0)&(y+yc<Cols)
        i(-x+xc,  y+yc) = value;
    end
end

%Fill the holes
i=imfill(i,'holes');


%Get rid of the boundaries used to fill regions that were on the border
i(1,:)=0;
i(:,1)=0;
i(Rows,:)=0;
i(:,Cols)=0;


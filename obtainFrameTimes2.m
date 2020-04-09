function obtainFrameTimes2(timeStampList, NSeries)

Frame_Times = [];

for s = 1:NSeries
   
   timeStampList = [];
   
   if s == 1
    first_time = hex2dec(timeStampList(1));
   end
   
   for f = 1:NFrames(s)
       FrameTimes(end+1) = hex2dec(timeStampList(f)) - first_time;
   end
    
end
function [hx, asc] = hexdump(filenm, n)
% hexdump(filenm, n)
% Print the first n bytes of a file in hex and ASCII.

fid = fopen(filenm, 'r');

if (fid<0), disp(['Error opening ',filenm]); return; end

hx = [];
asc = [];

nread = 0;

while (nread < n)
    width = 16;
    [A,count] = fread(fid, width, 'uchar');
    nread = nread + count;
    if (nread>n), count = count - (nread-n); A = A(1:count); end
    hexstring = repmat(' ',1, width*3);
    hexstring(1:3*count) = sprintf('%02x ',A);
    ascstring = repmat('.',1, count);
    idx = find(double(A)>=32);
    ascstring(idx) = char(A(idx));
    fprintf(1, '%s *%s*\n', hexstring, ascstring);
    hx = [hx, hexstring];
    asc = [asc, ascstring];
end

fclose(fid);

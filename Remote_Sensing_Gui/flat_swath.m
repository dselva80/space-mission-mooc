function [ A,minrow,maxrow,mincol,maxcol] = flat_swath( circ_rows,circ_cols,latg,r,height,width,B)
% OUTPUTS
% A where A is an image

% INPUTS
% B is an image

% set up some indexes in circle lon and circle lat
quart = round(length(circ_rows)/4);
half = round(length(circ_cols)/2);
tf = quart+half;

% this is the basecase (but useful to put it here)
minrow = min(circ_rows);
maxrow = max(circ_rows);
mincol = min(circ_cols);
maxcol = max(circ_cols);


if latg+r >= 90
    % at north pole
    % have a circle of rows and cols that cannot be fully displayed
    % add another flipped image to top of map so that the edges aren't a
    % problem. This can make the poles look a bit strange
    
    
    % have to split up the circle in half to add a displacement to only
    % half
    bottom_rows = circ_rows(quart+1:tf);
    top_rows = cat(1,circ_rows(1:quart),circ_rows(tf+1:end));
    top_rows=top_rows+(height-top_rows).*2;

    
    % map tries to display all cols here. it doesnt look great at large FOV
    
    minrow=min(bottom_rows);
    maxrow=max(top_rows);
    
    im = cat(1,B,fliplr(flip(B)));
    
    A = im(minrow:maxrow,mincol:maxcol);

 elseif latg-r<=-90
    % at south pole
    bottom_rows = circ_rows(quart+1:tf);
    top_rows = cat(1,circ_rows(1:quart),circ_rows(tf+1:end));
    bottom_rows = bottom_rows+height;
    top_rows = top_rows+height;

    
    
    minrow=min(bottom_rows);
    maxrow=max(top_rows);


    im = cat(1,flip(B,1),fliplr(B));

    A = im(minrow:maxrow,mincol:maxcol);
    
 elseif maxcol-mincol > width-10
    % at prime meridian must loop longitude

    beg_cols = zeros(size(circ_cols));
    end_cols = zeros(size(circ_cols));
    % split circular cols into two different arrays
    for q = 1:length(circ_cols)
        if circ_cols(q) < height
            beg_cols(q) = circ_cols(q);
        else
            end_cols(q) = circ_cols(q);
        end
    end
    beg_cols=beg_cols(beg_cols~=0)+width;
    end_cols=end_cols(end_cols~=0);

    maxcol = max(beg_cols);
    mincol = min(end_cols);

    im = cat(2,B,B);
    
    A = im(minrow:maxrow,mincol:maxcol);

else
    A = B(minrow:maxrow,mincol:maxcol);
end

end


function [ A,mincol, maxcol, minrow, maxrow ] = flat_swath( circ_rows,circ_cols,latg,r )
% OUTPUTS
% A where A is an image

% load iamge data
dat = load('topo.mat');
topolegend = dat.topolegend;
topo = dat.topo;

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
    % problem.
    
    
    % have to split up the circle in half to add a displacement to only
    % half
    bottom_rows = circ_rows(quart+1:tf);
    top_rows = cat(1,circ_rows(1:quart),circ_rows(tf+1:end));
    top_rows=top_rows+(180-top_rows).*2;

    
    % This is where I fudged it. It's not too noticeable at small (<50 FOV)
    % if cols are not fudged the map tries to display all cols and it looks
    % bad
    
    minrow=min(bottom_rows);
    maxrow=max(top_rows);
    
    mincol=150;
    maxcol=mincol+maxrow-minrow;


    im = cat(1,topo,fliplr(flip(topo)));

    A = im(minrow:maxrow,mincol:maxcol);


 elseif latg-r<=-90
    % at south pole
     bottom_rows = circ_rows(quart+1:tf);
    top_rows = cat(1,circ_rows(1:quart),circ_rows(tf+1:end));
    bottom_rows = bottom_rows+180;
    top_rows = top_rows+180;

    
    
    mincol=150;
    maxcol=200;
    minrow=min(bottom_rows);
    maxrow=max(top_rows);


    im = cat(1,flip(topo,1),fliplr(topo));

    A = im(minrow:maxrow,mincol:maxcol);
 elseif maxcol-mincol > 350
    % at prime meridian must loop longitude

    beg_cols = zeros(size(circ_cols));
    end_cols = zeros(size(circ_cols));
    % split circular cols into two different arrays
    for q = 1:length(circ_cols)
        if circ_cols(q) < 180
            beg_cols(q) = circ_cols(q);
        else
            end_cols(q) = circ_cols(q);
        end
    end
    beg_cols=beg_cols(beg_cols~=0)+360;
    end_cols=end_cols(end_cols~=0);

    maxcol = max(beg_cols);
    mincol = min(end_cols);

    im = cat(2,(topo),(topo));
    A = im(minrow:maxrow,mincol:maxcol);

else
    A = topo(minrow:maxrow,mincol:maxcol);
end


end

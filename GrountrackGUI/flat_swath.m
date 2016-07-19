function [ A, mincol, maxcol, minrow, maxrow ] = flat_swath( circ_rows,circ_cols,latg,r,height,width,indec )
% OUTPUTS
% A where A is an image

% load image data
dat = load('topo.mat');
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
    disp('north pole')
    % at north pole
    % have a circle of rows and cols that cannot be fully displayed
    % add another flipped image to top of map so that the edges aren't a
    % problem.
    
   
    
    % have to split up the circle in half to add a displacement to only
    % half
    
    % bottom rows is where x is neg
    bottom_rows = circ_rows(indec);
    circ_rows(indec)=400;
    top_rows = circ_rows(circ_rows~=400);
%     top_rows=height+(height-min(top_rows));
    maxrow = height+min(top_rows);
    minrow = min(circ_rows);
    
    

    if maxcol-mincol > width-20
        disp('both')
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
        beg_cols=beg_cols(beg_cols~=0)+width;
        end_cols=end_cols(end_cols~=0);

        maxcol = max(beg_cols);
        mincol = min(end_cols);
    end
    
    
    
    
%     if r>25
%         disp('r')
%         maxrow=height;
%     end
  

    im1 = cat(1,topo,rot90(topo,-2));
    im = cat(2,im1,im1);

    A = im(minrow:maxrow,mincol:maxcol);


 elseif latg-r<=-90
    % at south pole
    disp('south pole')
    top_rows = circ_rows(indec);
    circ_rows(indec) = 400;
    bottom_rows = circ_rows(circ_rows~=400);
    bottom_rows = bottom_rows+height;
    top_rows = height-top_rows;

    
    minrow=min(bottom_rows);
    maxrow=max(top_rows);

    

    im1 = cat(1,fliplr(flip(topo)),topo);
    im = cat(2,im1,im1);
    A = im(minrow:maxrow,mincol:maxcol);
    
 elseif maxcol-mincol > width-20
    % at prime meridian must loop longitude
    disp('prime')
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

    im = cat(2,topo,topo);
    A = im(minrow:maxrow,mincol:maxcol);

else
    A = topo(minrow:maxrow,mincol:maxcol);
end


end


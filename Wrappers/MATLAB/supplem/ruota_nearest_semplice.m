function [ out ] = ruota_nearest_semplice( img , gradi)
    % Rotates the image applying a nearest-neighbour interpolation
    % but computing for every pixel inside the bounding box
    [x, y]=size(img);
    gradi=(gradi/180)*pi;
    kernel=[cos(-gradi) sin(-gradi);-sin(-gradi) cos(-gradi)];
    c1=[x y]./2;
    d=ceil([ abs(x*cos(gradi)) + y*abs(sin(gradi)) abs(x*sin(gradi)) + abs(y*cos(gradi))]);
    out=ones(d);
    c2= round(d./2);
    for i=1:d(1) %x
        for j=1:d(2) %y
            p=([i j]-c2);
            p1=round(p*kernel)+c1;
            if(p1(1)>0 && p1(2)>0 && p1(1)<=x && p1(2)<=y)
                out(i,j)=img(p1(1),p1(2)); 
            end
        end
    end
end
function flood_fill(img,x,y,newcolor)

if(x>size(img,2)||y>size(img,1)||y<1||x<1)
    return
elseif(img(x,y)==1)
    img(x,y) = newcolor;
    %save('m.mat','img');
    imshow(img,[])
    drawnow
    flood_fill(img,x-1,y-1,newcolor)
    flood_fill(img,x,y-1,newcolor)
    flood_fill(img,x+1,y-1,newcolor)
    flood_fill(img,x-1,y,newcolor)
    flood_fill(img,x+1,y,newcolor)
    flood_fill(img,x-1,y+1,newcolor)
    flood_fill(img,x,y+1,newcolor)
    flood_fill(img,x+1,y+1,newcolor)
else
    
end

end
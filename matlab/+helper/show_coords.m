function show_coords(image, ROIs, CL, CR)

imagesc(image), colormap gray, hold on
for i = 1:size(CL,1)
    plot(CL(i,1),CL(i,2),'r*');
    %         plot(CR(i,1),CR(i,2),'r*');
    text(CR(i,1),CR(i,2),ROIs{i},'Color','red');
end
end
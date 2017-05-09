function showVideo4Test(sifImages)

for i=1:numel(sifImages)
    figure(1)
    imshow(sifImages{i}.image,[]);
    pause(0.01)
end
    
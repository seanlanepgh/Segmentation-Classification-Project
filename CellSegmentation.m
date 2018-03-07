%%
clc;	% Clear command window.
clear;	% Delete all variables.
close all;	% Close all figure windows except those created by imtool.
% imtool close all;	% Close all figure windows created by imtool.
workspace;	% Make sure the workspace panel is showing.
%%
srcFiles = dir('C:\Users\Sean Lane\Desktop\CellImages\*.tif');

for j = 1: length(srcFiles)
    filename = strcat('C:\Users\Sean Lane\Desktop\CellImages\',srcFiles(j).name);
     info = imfinfo(filename);
     disp(info)
   
    I = uint8(imread(filename));
     se = strel('disk',9);
    J = imsubtract(imadd(I,imtophat(I,se)),imbothat(I,se));
    % ... Do something with image A ...
   
%%
   

Equal = histeq(I,256);

%Contrast = imadjust(Equal);
%figure,imshow(Contrast);

    %Imed = medfilt2(I,[5 5], 'zeros');
    %figure,imshow(Imed);
%}
   
hy = fspecial('sobel');
hx = hy';
Iy = imfilter(double(J), hy, 'replicate');
Ix = imfilter(double(J), hx, 'replicate');
gradmag = sqrt(Ix.^2 + Iy.^2);
%figure ,imshow(gradmag,[]), title('Gradient magnitude (gradmag)')
L = watershed(gradmag);
Lrgb = label2rgb(L);
%figure, imshow(Lrgb), title('Watershed transform of gradient magnitude (Lrgb)')
se = strel('disk', 10);
Io = imopen(J, se);
%figure
%imshow(Io), title('Opening (Io)')
  Ie = imerode(J, se);
Iobr = imreconstruct(Ie,J);
%figure
%imshow(Iobr), title('Opening-by-reconstruction (Iobr)') 
Ioc = imclose(Io, se);
%figure
%imshow(Ioc), title('Opening-closing (Ioc)')
Iobrd = imdilate(Iobr, se);
Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);
%figure
%imshow(Iobrcbr), title('Opening-closing by reconstruction (Iobrcbr)')
fgm = imregionalmax(Iobrcbr);
%figure
%imshow(fgm), title('Regional maxima of opening-closing by reconstruction (fgm)')
I2 = I;
I2(fgm) = 255;
%figure
%imshow(I2), title('Regional maxima superimposed on original image (I2)')
se2 = strel('disk',5);
fgm2 = imclose(fgm, se2);
fgm3 = imerode(fgm2, se2);
fgm4 = bwareaopen(fgm3, 20);
I3 = I;
I3(fgm4) = 255;
%figure
%imshow(I3)
%title('Modified regional maxima superimposed on original image (fgm4)')
bw = imbinarize(Iobrcbr);
%figure
%imshow(bw), title('Thresholded opening-closing by reconstruction (bw)')
%D = -bwdist(~bw,'cityblock');
D=bwdist(bw);
%figure, imshow(D,[])
DL = watershed(D,8);
bgm = DL == 0;
%figure
%imshow(bgm), title('Watershed ridge lines (bgm)')
bw2 =bw;
bw2(DL == 0) =0;
%figure,imshow(bw2);
mask = imextendedmin(D,2);
%imshowpair(bw,mask,'blend')
D2 = imimposemin(D,mask);
Ld2 = watershed(D2);
bw3 = bw;
bw3(Ld2 == 0) = 0;
%figure ,imshow(bw3)
gradmag2 = imimposemin(gradmag, bgm | fgm4);
L = watershed(gradmag2);
I4 = I;
I4(imdilate(L == 0, ones(3, 3)) | bgm | fgm4) = 255;
%figure
%imshow(I4)
%title('Markers and object boundaries superimposed on original image (I4)')
Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
%figure
%imshow(Lrgb)
%title('Colored watershed label matrix (Lrgb)')
%figure
%imshow(I)
%hold on
%himage = imshow(Lrgb);
%himage.AlphaData = 0.3;
%title('Lrgb superimposed transparently on original image')
 bw_perim = bwperim(L);
                    overlay = imoverlay(I,bw_perim,[1 .3 .3]);
imwrite(overlay,(['\Users\Sean Lane\Desktop\Watershed overlay\',srcFiles(j).name]));
%}
end
   
  

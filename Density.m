clear; close all; clc;
%Disable 'too big' warning
WarnID = 'images:initSize:adjustingMag';
warning('off',WarnID);

%import mammogram
Mam = imread('mdb005.pgm');
MamComp = imread('mdb006.pgm');

Mam = imsharpen(Mam);
figure;
imshow(Mam);
title('1-Right')

MamComp = imsharpen(MamComp);
figure;
imshow(MamComp);
title('2-Left')

%Use below functions to isolate breast and estimate tissue density
%Right
Mam = Isolate(Mam);

disp('Percentage density = ' )
MamTissue = Tissue(Mam);

figure;
imshow(Mam);
title('3-Right')

%Left
MamComp = Isolate(MamComp);

disp('Percentage density = ' ) 
MamCompTissue = Tissue(MamComp);

figure;
imshow(MamComp);
title('4-Left')

%if there is a tissue density difference between the two breast ensities, flag a
%problem
if(MamCompTissue - MamTissue > 0.2 || MamTissue - MamCompTissue > 0.2)
    BIRADs = 2;
    
    MamTumourR = FindTumour(Mam);
    
    figure;
    imshow(MamTumourR);
    title('5-Right')

    MamTumourL = FindTumour(MamComp);
    
    figure;
    imshow(MamTumourL);
    title('6-Left')
%Tumor finding through comparrison of breasts and edge
    FinalTumour = Symetrical(MamTumourR, MamTumourL);
else
    BIRADs = 1;
end

%Label shapes to compare them
Perimeter = regionprops(FinalTumour, 'perimeter');
Perimeter = struct2array(Perimeter);
Area = regionprops(FinalTumour, 'area');
Area = struct2array(Area);
Pos = regionprops(FinalTumour, 'Centroid');
Pos = struct2array(Pos);

%if shape is detected, BIRADS = 3
if(isempty(Perimeter) == false)
    BIRADs = 3;
end

MamComp = flip(MamComp, 2);

for L = 1:length(Area)
% Compute roundness and designare tumor type. Give TNM value and position
    if(Pos(1,L) > 1024)
        Pos(1,L) = Pos(1,L) - 1024;
    end
    Round(1,L) = 4*pi*Area(1,L)/Perimeter(1,L)^2;
    disp('Tumor ID =');
    disp(L);
    disp('X Co-Ordinate =');
    disp(Pos(1,L));
    disp('Y Co-Ordinate =');
    disp(Pos(1,(L*2)));
    disp('BIRADS value =');
    disp(BIRADs);
    disp('Tumour roundness =');
    disp(Round(1,L));
    disp('Tumour area =');
    disp(Area(1,L)); 
    disp('-------------------------------------------');

end

MamCompBW = FinalTumour(:, 1:1024, :);
MamBW = FinalTumour(:, 1025:end, :);


MamCompFin = imfuse(MamComp, MamCompBW);
MamFin = imfuse(Mam, MamBW);

figure;
imshow(MamCompFin);
title('Left_Final');
figure;
imshow(MamFin);
title('Right_Final');

%%
%%============================Isolating breast tissue===============================
%Get rid of external noise and objects e.g. identifying tag
%Convert to binary and areaopen
function Mam = Isolate(Mam)

MamBin = imbinarize(Mam,0.1);
for X = 1:1024
    MamBin(1,X) = 0;
end
MamBin = bwareaopen(MamBin,100000);

RightX = 1024;
LeftX = 1;
%Check image dimensions
for LeftX = 2:1024
    RightX = RightX - 1;
    if(MamBin(2,LeftX) > 0)
        break
    else
        if(MamBin(2,RightX) > 0)
            break
        end
        
    end
end
RightX = RightX-1;
LeftX = LeftX+1;

%Any pixels that are converted to 0 in the binary image are converted to 0
%in the original mammogram
for Y = 1:1024
    for X = 1:1024
        if(MamBin(X,Y) == 0)
            Mam(X,Y) = 0;
        end
    end
end


%Identify orientation of pectoral muscle
if(Mam(2, (LeftX)) > 0 && LeftX < 500)
    disp('Left Peck')
    for X = 2:1024
        for Y = LeftX:1024
            %Run through top-left of screen changing all dark index values
            %to 0 until the lighter area of breast is reached
            if(Mam(X,Y) < 190)
            break
            end
            Mam(X,Y) = 0;
        end
    end
else  
    if(Mam(2, (RightX)) > 0 && RightX > 700)
        disp('Right Peck')
        for X = 2:1024
            for Y = RightX:-1:1
                if(Mam(X,Y) < 190)
                break
                end
                Mam(X,Y) = 0;
            end
        end
    else
        disp('No pectoral muscle - place differently')
    end
end

end


%%
%%===========================Density Calculation============================
%Calculating the percentage of the breast with a darkness above set threshold
function DensityPercent = Tissue(Mam)

FullCount = 0;
DenseCount = 0;
for X = 1:1024
    for Y = 1:1024
        if(Mam(X,Y) ~= 0)
           FullCount = FullCount + 1;
           if(Mam(X,Y) > 170)
               DenseCount = DenseCount + 1;
           end
        end
    end
end


DensityPercent = (DenseCount/FullCount)*100;

disp(DensityPercent)
disp('Tissue is:')
if(DensityPercent <= 10)
    Dense = 1;
    disp('A: Fatty' )
else
    if(DensityPercent > 10 && DensityPercent <= 25)
        Dense = 2;
        disp('B: Scattered Fibroglandular' )
    else
        if(DensityPercent > 25 && DensityPercent <= 45)
            Dense = 3;
            disp('C: Heterogeneously Dense')
        else
            if(DensityPercent > 45)
                Dense = 4;
                disp('D: Dense')
            end
        end
    end
end

end
%%
%%============================Separating Tumour===========================
function MamEdge = FindTumour(Mam)
%Bluring image with gaussian filtering
Mam = imgaussfilt(Mam,5);


%Finding edge using canny
Mamedge = edge(Mam, 'canny', 0.1);
%apply morphology to thicken lines and close gaps
Nhood = strel('disk', 3);
Mamedge = imdilate(Mamedge, Nhood);
%Thin all lines to their skeleton without decreasing length; similar to
%eroding but maintains length of lines - better effect in this case
Mamedge = bwmorph(Mamedge, 'thin', inf);
Mamedge = bwmorph(Mamedge, 'diag');



%Connected components to identify all lines and shapes in edge image
MamedgeCC = bwconncomp(Mamedge);

%Closing of shapes through dilating in a loop

BWblank = false(MamedgeCC.ImageSize);
stats = regionprops(MamedgeCC,'ConvexImage','EulerNumber');
%calculates new image using euclidian distance transformation of Img to
%close major gaps
for Img = find([stats.EulerNumber]>0)
    DistImg = bwdist(~stats(Img).ConvexImage);
    MaxClosed = ceil(max(DistImg(:)));
    BWSlice = BWblank;
    BWSlice(MamedgeCC.PixelIdxList{Img}) = true;
    if isinf(MaxClosed)
        continue; 
    end
    for dilSz = 2:MaxClosed
        BWNew = imdilate(BWSlice,ones(dilSz));
        StatsNew = regionprops(BWNew,'EulerNumber');
        if StatsNew.EulerNumber<=0
            BWNew = imerode(imfill(BWNew,'holes'),ones(dilSz));
            MamedgeCC.PixelIdxList{Img} = find(BWNew);
            %if number of objects-number of holes is <= 0, erode that
            %object and fill any holes
        end
    end
end
MamEdge = imfill(labelmatrix(MamedgeCC),'holes');

stats = regionprops(MamedgeCC,'ConvexImage','EulerNumber','BoundingBox');
for Img = find([stats.EulerNumber]>0)
    MaxClosed = ceil(max(DistImg(:)));
    BWSlice = BWblank;
    BWSlice(MamedgeCC.PixelIdxList{Img}) = true;
    DistImg = bwdist(~BWSlice);
    if ~any(DistImg(:)>1)
        BWNew = BWSlice;
        BoundBox = ceil(stats(Img).BoundingBox);
        BWNew((1:BoundBox(4))+BoundBox(2)-1,(1:BoundBox(3))+BoundBox(1)-1) = stats(Img).ConvexImage;
        MamedgeCC.PixelIdxList{Img} = find(BWNew);
    end
end
Mam = imfill(labelmatrix(MamedgeCC),'holes');

end

% % ========================Finding unsymetrical areas===================

%%
function Diff = Symetrical(Mam, MamComp)

%flip one image
MamComp = flip(MamComp,2);
CCR = Mam;
CCL = MamComp;

%%Right side
%Compute area for each white shape
AreaR = regionprops(CCR,'area');
AreaR = struct2array(AreaR);
%Compute perimeter for each white shape
PerimeterR = regionprops(CCR,'Perimeter');
PerimeterR = struct2array(PerimeterR);

RoundR(1,1) = 256;
RoundR(1,2) = 256;
for R = 3:length(AreaR)
%Compute roundness matrix for each object and eliminate any below a given
%threshold
    
    
    RoundR(1,R) = 4*pi*AreaR(1,R)/PerimeterR(1,R)^2;
    if(RoundR(1,R) > 0.4 && PerimeterR(1,R) > 20 && AreaR(1,R) < 10000)
        RoundR(1,R) = R;
    else 
        RoundR(1,R) = 256;
    end
   

end
CCRNew = ismember(CCR, RoundR);

figure;
imshow(CCRNew);
title('Right')

%again for left 

%Compute area for each white shape
AreaL = regionprops(CCL,'area');
AreaL = struct2array(AreaL);
%Compute perimeter for each white shape
PerimeterL = regionprops(CCL,'Perimeter');
PerimeterL = struct2array(PerimeterL);

RoundL(1,1) = 256;
RoundL(1,2) = 256;

for L = 3:length(AreaL)
%compute roundness matrix for each object and eliminate any below a given
%threshold
    RoundL(1,L) = 4*pi*AreaL(1,L)/PerimeterL(1,L)^2;
    if(RoundL(1,L) > 0.4 && PerimeterL(1,L) > 20 && AreaL(1,L) < 10000)
        RoundL(1,L) = L;
    else 
        RoundL(1,L) = 256;
    end
end

CCLNew = ismember(CCL, RoundL);

figure;
imshow(CCLNew);
title('Left');

%Label shapes to compare them
PerimeterL = regionprops(CCLNew, 'perimeter');
PerimeterL = struct2array(PerimeterL);
PerimeterR = regionprops(CCRNew, 'perimeter');
PerimeterR = struct2array(PerimeterR);
AreaL = regionprops(CCLNew, 'area');
AreaL = struct2array(AreaL);
AreaR = regionprops(CCRNew, 'area');
AreaR = struct2array(AreaR);

%Comparing shapes through their perimeters and areas
%if Area L is within +-20 of AreaR and Perimeter R is within +- 20 of
%perimeter L, delete both
for Y = 1:length(AreaL)
    for X = 1:length(AreaR)
        %if the anomoly is identical on both sides, it is likely to be
        %normal
        if(sqrt((AreaL(Y) - AreaR(X))^2) <= 20 && sqrt((PerimeterL(Y) - PerimeterR(X))^2) <= 20)
            AreaL(Y) = 0;
            AreaR(X) = 0;
            CCLNew = ismember(CCLNew,AreaL);
            CCRNew = ismember(CCRNew, AreaR);

        end
    end
end

Diff = [CCLNew, CCRNew];





end



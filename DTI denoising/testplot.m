clear all
close all
addpath example_dwi
addpath NIfTI_20140122
addpath fanDTasia
nii_file = load_untouch_nii('nifti_dt.nii');
data = nii_file.img;
clear nii_file;
% for k = 1:size(data,3)
k = 28;
slice = data(:,:,k,:,:);
slice2 = reshape(slice,112,112,6);

tensor = zeros(3,3,112,112);

for i = 1:112
    for j = 1:112
        tensor(1,1,i,j) = slice2(i,j,1);
        tensor(2,1,i,j) = slice2(i,j,2);
        tensor(1,2,i,j) = slice2(i,j,2);
        tensor(2,2,i,j) = slice2(i,j,3);
        tensor(3,1,i,j) = slice2(i,j,4);
        tensor(1,3,i,j) = slice2(i,j,4);
        tensor(3,2,i,j) = slice2(i,j,5);
        tensor(2,3,i,j) = slice2(i,j,5);
        tensor(3,3,i,j) = slice2(i,j,6);
    end
end

tensor = tensor./max(max(max(max(tensor))));
tensor = tensor*5;
tic
close all
plotDTI(tensor(:,:,28:87,34:73))
toc
% pause(5)
% end
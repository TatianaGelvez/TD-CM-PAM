function Sgnl = A_direct_conv2D(par)
circKernel_flip = flip(flip(par.circKernel,1),2);
Sgnl_Conv2D = 0;
my_image        = reshape(par.BB, par.NPixlsX, par.NPixlsZ, par.NSampls);
my_image        = permute(my_image, [2 1 3]);
my_image        = cat(3, my_image, my_image);   % (Z × X × 2*Ns)

if (par.AcqType == "ideal")
for con = 1:par.NPixlsZ
    imageSlice  = squeeze(my_image(con,:,:)); 
    kernelSlice = squeeze(circKernel_flip(end - con + 1,:,:));
    Sgnl_Conv2D = Sgnl_Conv2D + conv2(imageSlice, kernelSlice, "full"); 
end
else
for con = 1:par.NPixlsZ
    imageSlice  = squeeze(my_image(con,:,:)); 
    kernelSlice = squeeze(circKernel_flip(end - con + 1,:,:));
    Sgnl_Conv2D = Sgnl_Conv2D + conv2(imageSlice, kernelSlice, "full")./(sqrt(par.NPixls*par.NSensor)); 
end
end

%%Extract relevant output
CenterResultX = round(size(Sgnl_Conv2D,1)/2);
Sgnl    = Sgnl_Conv2D(CenterResultX-par.centerSensor+1:CenterResultX+par.centerSensor, par.NSampls+1:2*par.NSampls);
end

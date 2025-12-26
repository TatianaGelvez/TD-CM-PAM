function BB = A_transp_conv2D(par)

Y_padded       = zeros(par.NPixlsX + size(par.circKernel,2) - 1 ,   par.NSampls*3-1);
CenterResultX  = round(size(Y_padded,1) / 2);
Y_padded(CenterResultX-par.centerSensor+1:CenterResultX+par.centerSensor, par.NSampls+1:2*par.NSampls) = par.Sgnl;

my_image_est   = zeros(par.NPixlsZ, par.NPixlsX, 2*par.NSampls);  % (Z × X × 2*Ns)

for con = 1:par.NPixlsZ
    kernelSlice_adj = squeeze(par.circKernel_adj(con,:,:));
    my_image_est(con,:,:) = conv2(Y_padded, kernelSlice_adj, 'valid')./(sqrt(par.NPixls*par.NSensor));
end
BB = my_image_est(:,:,1:par.NSampls);
BB = BB + my_image_est(:,:,par.NSampls + 1:2*par.NSampls);
BB = permute(BB, [2 1 3]);
BB = reshape(BB, par.NPixls, par.NSampls);
end

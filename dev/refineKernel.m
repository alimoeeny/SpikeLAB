function [ newKernel ] = refineKernel( oldKernel, intendedLuminances, measuredLuminances )
% refine Kernel gets the old gamma lookup kernel (vector pf 256 uint8) and
% the measured Luminances (vector of 3 or more elements can be doubles) and
% generates a new Kernel that if used instead of the old Kernel generates a
% linear output
% assumption here is that there are 3 gamma functions:
% 1- the native gamma function of the device
% 2- the old gamma function applied to do the measuments
% 3- the measured output which is also a gamma function
% old(native(x)) = measured(x)
% x ^ (oldGamma * nativeGamma) = x ^ (measuredGamma)
% nativeGamma = measuredGamma / oldGamma;
% newKernel = [1:1:256] .^ (1/nativeGamma)

assert(sum(size(oldKernel)==[1 256])==2);
assert(sum(size(measuredLuminances)==size(intendedLuminances))==2);

[cfun,gof,output] = fit(intendedLuminances', measuredLuminances', fittype('a*(x^b)'));

measuredGamma = cfun.b;
mscalefactor = cfun.a;

[cfun,gof,output] = fit([1:1:256]', oldKernel', fittype('a*(x^b)'));
oldGamma = cfun.b
oscalefactor = cfun.a;

newKernel = [1:1:256] .^ (1/(measuredGamma/oldGamma));





end


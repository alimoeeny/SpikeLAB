clear;
w = 2000;
h = 2000;
panx = -2250; %-3300;
pany = -680;% -1300;
constant = 2;
zoom = 3000; %5965;

I = zeros(w,h);

parfor a = 1:w
    for b = 1:h
        c = ((a+panx-(w/2))+(b+pany-(h/2))*i) / zoom;
        v = 0;
        z = 0i+0;
        for t = 1:100000
            z = (z^2) + c;
        end
        v = abs(z^200);
        I(a,b) = log(v);
    end
end

cd ~/Dropbox//Projects/SpikeLAB/
imwrite(I, 'fractal.bmp', 'bmp');

figure(110), imagesc(I);
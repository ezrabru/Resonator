function rgbTriplet = getRGBtripletFromWavelength(wavelength_nm)

col800 = [0.4 0.0 0.0];
col700 = [0.7 0.0 0.0];
col650 = [1.0 0.0 0.0];
col630 = [1.0 0.1 0.0];
col561 = [0.9 0.9 0.0];
col515 = [0.6 0.9 0.0];
col488 = [0.0 0.9 0.9];
col475 = [0.0 0.7 0.9];
col450 = [0.0 0.3 0.9];
col405 = [0.8 0.0 1.0];

wavelength = [405, 450, 475, 488, 515, 561, 630, 650, 700, 800];
rgb = [col405;
       col450;
       col475;
       col488;
       col515;
       col561;
       col630;
       col650;
       col700;
       col800];
r = rgb(:,1);
g = rgb(:,2);
b = rgb(:,3);

x = wavelength;
xq = 400:800;
method = 'makima';
rq = interp1(x,r,xq,method); rq(rq<0) = 0; rq(rq>1) = 1; rq = rq';
gq = interp1(x,g,xq,method); gq(gq<0) = 0; gq(gq>1) = 1; gq = gq';
bq = interp1(x,b,xq,method); bq(bq<0) = 0; bq(bq>1) = 1; bq = bq';

[~,idx_wavelength] = min(abs(xq - wavelength_nm));
rgbTriplet = [rq(idx_wavelength) gq(idx_wavelength) bq(idx_wavelength)];

% plot(x,r,'o',xq,rq,':.'); hold on
% plot(x,g,'o',xq,gq,':.');
% plot(x,b,'o',xq,bq,':.');

% spectrumImage = cat(3,rq',gq',bq');
% imagesc(spectrumImage,'XData',xq)

end


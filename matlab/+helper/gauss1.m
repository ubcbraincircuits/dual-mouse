function Y = gauss1(X,sigma,sz)
x = linspace(-sz / 2, sz / 2, sz);
gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
gaussFilter = gaussFilter / sum (gaussFilter); % normalize

Y = conv(X,gaussFilter,'same');
end


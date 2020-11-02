function roundans = roundk(input, k)
mult = 10^k;
roundans = round((input*mult))/mult;

end
function normal = checkNormal(data)
    normal = ~kstest(zscore(data));
end
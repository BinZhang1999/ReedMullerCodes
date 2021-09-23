function hadamardSetting = HadamardSetting(m)
hadamardSetting.m = m;
hadamardSetting.K = m + 1;
hadamardSetting.N = 2^m;
hadamardSetting.hadamardMatrix = hadamard(hadamardSetting.N);

hadamardSetting.G = zeros(m+1, hadamardSetting.N);
hadamardSetting.G(m:-1:1,:) = (myDec2Bin((0:(2^m-1))',m))';
hadamardSetting.G(m+1,:) = ones(1,hadamardSetting.N);
end
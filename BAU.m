clc;
clear;

load("2017-new-d.mat")
load("carbon.mat")
load("S1.mat")

nSec = 70;
nReg = 160;
N = nSec * nReg;

brandRows = [29 30 31 32 33];
sheetNames = {'cc-off','cc-hm','cc-zara','cc-gap','cc-uniqlo'};
salesRaw = {soff, shm, szara, sgap, suniqlo};

years = [2020 2025 2030 2035 2040 2045 2050];

sfList = {sf2020, sf2025, sf2030, sf2035, sf2040, sf2045, sf2050};
pfList = {pf2020, pf2025, pf2030, pf2035, pf2040, pf2045, pf2050};

salesVec = cell(1, length(brandRows));

for b = 1:length(brandRows)
    M = zeros(nSec, nReg);
    M(brandRows(b), :) = salesRaw{b};
    salesVec{b} = reshape(M, N, 1);
end

F = zeros(N, nReg);

for i = 1:nReg
    F(:, i) = sum(f3(:, (i-1)*3+1 : i*3), 2);
end

Ybase = cell(1, length(brandRows));

for b = 1:length(brandRows)
    Y = zeros(N, nReg);

    for r = 1:nReg
        rowID = brandRows(b) + (r-1) * nSec;
        Y(rowID, :) = F(rowID, :);
    end

    Ybase{b} = Y;
end

cc = reshape(carbon, N, 1);
cff = cc ./ x3;
cff(~isfinite(cff)) = 0;

X1 = repmat(x3', N, 1);
A = X ./ X1;
A(~isfinite(A)) = 0;

L2017 = inv(eye(N) - A);
L2017(~isfinite(L2017)) = 0;

cBase = cell(1, length(brandRows));
fBase = cell(1, length(brandRows));

for b = 1:length(brandRows)
    S = repmat(salesVec{b}, 1, nReg);

    cBase{b} = Ybase{b} ./ S;
    cBase{b}(~isfinite(cBase{b})) = 0;

    fBase{b} = S ./ p17;
    fBase{b}(~isfinite(fBase{b})) = 0;
end

for y = 1:length(years)

    yearNow = years(y);
    sf = sfList{y};
    pf = pfList{y};

    pYear = p17 .* pf;

    outFile = sprintf("s1c%02d-0305.xlsx", yearNow - 2000);

    for b = 1:length(brandRows)

        fYear = fBase{b} .* sf;

        Yyear = cBase{b} .* fYear .* pYear;
        Yyear(~isfinite(Yyear)) = 0;

        CFF = cff .* (L2017 * Yyear);

        CFC = squeeze(sum(reshape(CFF, nSec, nReg, nReg), 1));

        xlswrite(outFile, CFC, sheetNames{b}, "C3");

    end

end
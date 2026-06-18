clc;
clear;

load("2017-new-d.mat")
load("carbon.mat")
load("S1.mat")

nSec = 70;
nReg = 160;
N = nSec * nReg;

brandRows = [29 30 31 32 33];
sheetNames = {'cc-Other fast fashion','cc-H&M','cc-Inditex','cc-Gap','cc-Fast Retailing'};
salesRaw = {sOther fast fashion, sH&M, sInditex, sGap, sFast Retailing};

years = [2020 2025 2030 2035 2040 2045 2050];

sfList = {sf2020, sf2025, sf2030, sf2035, sf2040, sf2045, sf2050};
pfList = {pf2020, pf2025, pf2030, pf2035, pf2040, pf2045, pf2050};

tpList = {tp2020c, tp2025c, tp2030c, tp2035c, tp2040c, tp2045c, tp2050c};
fcList = {f2020c, f2025c, f2030c, f2035c, f2040c, f2045c, f2050c};

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

    Y = zeros(N, nReg);
    idx = brandRows(b):nSec:N;
    Y(idx, :) = F(idx, :);

    S = repmat(salesVec{b}, 1, nReg);

    cBase{b} = Y ./ S;
    cBase{b}(~isfinite(cBase{b})) = 0;

    fBase{b} = S ./ p17;
    fBase{b}(~isfinite(fBase{b})) = 0;

end

for y = 1:length(years)

    yearNow = years(y);

    sf = sfList{y};
    pf = pfList{y};
    tp = tpList{y};
    fc = fcList{y};

    tpMat = repmat(tp', nSec, 1);
    fcMat = repmat(fc, 1, nReg);

    cfYear = reshape(tpMat .* fcMat, N, 1) .* cff;
    cfYear(~isfinite(cfYear)) = 0;

    pYear = p17 .* pf;

    outFile = sprintf("s3c%02d-0305.xlsx", yearNow - 2000);

    for b = 1:length(brandRows)

        fYear = fBase{b} .* sf;

        Yyear = cBase{b} .* fYear .* pYear;
        Yyear(~isfinite(Yyear)) = 0;

        CFF = cfYear .* (L2017 * Yyear);

        CFC = squeeze(sum(reshape(CFF, nSec, nReg, nReg), 1));

        xlswrite(outFile, CFC, sheetNames{b}, "C3");

    end

end
clc;
clear;

load("2017-new-d.mat")
load("carbon.mat")
load("S5.mat")
load("mfuture.mat")

nSec = 70;
nReg = 160;
N = nSec * nReg;

brandNames = {'Other fast fashion','H&M','Inditex','gap','Fast Retailing'};
brandRows = [29 30 31 32 33];
sheetNames = {'cc-Other fast fashion','cc-H&M','cc-Inditex','cc-gap','cc-Fast Retailing'};

salesRaw = {sOther fast fashion, sH&M, sInditex, sgap, sFast Retailing};
mRaw = {mOther fast fashion, mH&M, mInditex, mgap, mFast Retailing};

years = [2020 2025 2030 2035 2040 2045 2050];
yearTags = [20 25 30 35 40 45 50];

sfList = {sf2020, sf2025, sf2030, sf2035, sf2040, sf2045, sf2050};
pfList = {pf2020, pf2025, pf2030, pf2035, pf2040, pf2045, pf2050};

salesVec = cell(1, length(brandRows));
mVec = cell(1, length(brandRows));

for b = 1:length(brandRows)
    M = zeros(nSec, nReg);
    M(brandRows(b), :) = salesRaw{b}(:)';
    salesVec{b} = reshape(M, N, 1);

    M = zeros(nSec, nReg);
    M(brandRows(b), :) = mRaw{b}(:)';
    mVec{b} = reshape(M, N, 1);
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
mBase = cell(1, length(brandRows));

for b = 1:length(brandRows)

    Y = zeros(N, nReg);
    idx = brandRows(b):nSec:N;
    Y(idx, :) = F(idx, :);

    S = repmat(salesVec{b}, 1, nReg);

    cBase{b} = Y ./ S;
    cBase{b}(~isfinite(cBase{b})) = 0;

    fBase{b} = S ./ p17;
    fBase{b}(~isfinite(fBase{b})) = 0;

    mBase{b} = repmat(mVec{b}, 1, nReg);

end

for y = 1:length(years)

    yearNow = years(y);
    yearTag = yearTags(y);

    sf = sfList{y};
    pf = pfList{y};

    pYear = p17 .* pf;

    outFile = sprintf("s4c%02d-0305.xlsx", yearNow - 2000);

    for b = 1:length(brandRows)

        varName = sprintf("m%s%d", brandNames{b}, yearTag);
        mFutureRaw = eval(varName);

        M = zeros(nSec, nReg);
        M(brandRows(b), :) = mFutureRaw(:)';
        mFutureVec = reshape(M, N, 1);
        mFuture = repmat(mFutureVec, 1, nReg);

        rYear = mFuture ./ mBase{b};
        rYear(~isfinite(rYear)) = 0;

        fYear = fBase{b} .* sf;

        Yyear = cBase{b} .* fYear .* pYear .* rYear;
        Yyear(~isfinite(Yyear)) = 0;

        CFF = cff .* (L2017 * Yyear);

        CFC = squeeze(sum(reshape(CFF, nSec, nReg, nReg), 1));

        xlswrite(outFile, CFC, sheetNames{b}, "C3");

    end

end
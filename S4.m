clc;
clear;

load("2017-new-d.mat")
load("carbon.mat")
load("S5.mat")
load("mfuture.mat")

nSec = 70;
nReg = 160;
N = nSec * nReg;

brandNames = {'off','hm','zara','gap','uniqlo'};
brandRows = [29 30 31 32 33];
sheetNames = {'cc-off','cc-hm','cc-zara','cc-gap','cc-uniqlo'};

salesRaw = {soff, shm, szara, sgap, suniqlo};
mRaw = {moff, mhm, mzara, mgap, muniqlo};

years = [2020 2025 2030 2035 2040 2045 2050];
yearTags = [20 25 30 35 40 45 50];

sfList = {sf2020, sf2025, sf2030, sf2035, sf2040, sf2045, sf2050};
pfList = {pf2020, pf2025, pf2030, pf2035, pf2040, pf2045, pf2050};

tpList = {tp2020c, tp2025c, tp2030c, tp2035c, tp2040c, tp2045c, tp2050c};
fcList = {f2020c, f2025c, f2030c, f2035c, f2040c, f2045c, f2050c};

salesVec = cell(1, length(brandRows));
mBaseVec = cell(1, length(brandRows));

for b = 1:length(brandRows)
    salesVec{b} = putRow(salesRaw{b}, brandRows(b), nSec, nReg);
    mBaseVec{b} = putRow(mRaw{b}, brandRows(b), nSec, nReg);
end

F = squeeze(sum(reshape(f3, N, 3, nReg), 2));

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

    mBase{b} = repmat(mBaseVec{b}, 1, nReg);

end

for y = 1:length(years)

    yearNow = years(y);
    yearTag = yearTags(y);

    sf = sfList{y};
    pf = pfList{y};
    tp = tpList{y};
    fc = fcList{y};

    tpMat = repmat(tp(:)', nSec, 1);
    fcMat = repmat(fc(:), 1, nReg);

    cfYear = reshape(tpMat .* fcMat, N, 1) .* cff;
    cfYear(~isfinite(cfYear)) = 0;

    pYear = p17 .* pf;

    outFile = sprintf("s5c%02d-0305.xlsx", yearNow - 2000);

    for b = 1:length(brandRows)

        mFutureName = sprintf("m%s%d", brandNames{b}, yearTag);
        mFutureRaw = eval(mFutureName);
        mFutureVec = putRow(mFutureRaw, brandRows(b), nSec, nReg);
        mFuture = repmat(mFutureVec, 1, nReg);

        rYear = mFuture ./ mBase{b};
        rYear(~isfinite(rYear)) = 0;

        fYear = fBase{b} .* sf;

        Yyear = cBase{b} .* fYear .* pYear .* rYear;
        Yyear(~isfinite(Yyear)) = 0;

        CFF = cfYear .* (L2017 * Yyear);

        CFC = squeeze(sum(reshape(CFF, nSec, nReg, nReg), 1));

        xlswrite(outFile, CFC, sheetNames{b}, "C3");

    end

end

function vec = putRow(v, rowID, nSec, nReg)
    M = zeros(nSec, nReg);
    M(rowID, :) = v(:)';
    vec = reshape(M, nSec * nReg, 1);
end
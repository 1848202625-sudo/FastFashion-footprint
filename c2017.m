clc
clear

load('2017-new-d.mat')
load('carbon.mat')

nSec = 70;
nReg = 160;
N = nSec * nReg;

outFile = "0927cd-17-160.xlsx";
rFile = "r2017.xlsx";

groupNames = {'all','traditional apperal','Other fast fashion','H&M','Inditex','gap','FastRetailing'};
prodSheets = {'production','traditional apperal','Other fast fashion','H&M','Inditex','gap','FastRetailing'};
ccSheets = {'cc-all','cc-traditional apperal','cc-Other fast fashion','cc-H&M','cc-Inditex','cc-gap','cc-FastRetailing'};
rSheets = {'c-traditional apperal','c-Other fast fashion','c-H&M','c-Inditex','c-gap','c-FastRetailing'};

groupRows = {
    28:33, ...
    28, ...
    29, ...
    30, ...
    31, ...
    32, ...
    33
};

cc = reshape(carbon, N, 1);
cff = cc ./ x3;
cff(~isfinite(cff)) = 0;

X1 = repmat(x3', N, 1);
A = X ./ X1;
A(~isfinite(A)) = 0;

L2017 = inv(eye(N) - A);
L2017(~isfinite(L2017)) = 0;

for g = 1:length(groupRows)

    Yg = zeros(N, 480);

    for r = 1:nReg
        rows = groupRows{g} + (r-1) * nSec;
        Yg(rows, :) = f3(rows, :);
    end

    CFF = cff .* (L2017 * Yg);

    cfp = sum(CFF, 2);
    cfpp = reshape(cfp, nSec, nReg);

    xlswrite(outFile, cfpp, prodSheets{g}, "C3");

    cfppp = squeeze(sum(reshape(CFF, nSec, nReg, 480), 1));

    cfpppp = squeeze(sum(reshape(cfppp, nReg, 3, nReg), 2));

    xlswrite(outFile, cfpppp, ccSheets{g}, "C3");

    if g >= 2
        xlswrite(rFile, cfppp, rSheets{g-1}, "C3");
    end

end
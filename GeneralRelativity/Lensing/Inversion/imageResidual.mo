within GeneralRelativity.Lensing.Inversion;
function imageResidual "How well a simulated image matches an observed one: sum of squares, chi-squared, explained fraction and correlation"
  input Real observed[:, :];
  input Real simulated[size(observed, 1), size(observed, 2)];
  input Real noiseSigma = 1.0 "Pixel noise standard deviation for chi-squared";
  output Real sumSquares "Sum over pixels of (observed - simulated)^2";
  output Real chi2 "sumSquares / noiseSigma^2; about the number of pixels when the model is right";
  output Real explained "1 - sumSquares / sum(observed^2): the fraction of the image's energy the simulation reproduces";
  output Real correlation "Pearson correlation of the two images, 1 for a perfect match up to scale and offset";
protected
  Integer n = size(observed, 1) * size(observed, 2);
  Real meanO = sum(observed) / n;
  Real meanM = sum(simulated) / n;
  Real sOO = 0.0;
  Real sMM = 0.0;
  Real sOM = 0.0;
  Real sO2 = 0.0;
algorithm
  sumSquares := 0.0;
  for i in 1:size(observed, 1) loop
    for j in 1:size(observed, 2) loop
      sumSquares := sumSquares + (observed[i, j] - simulated[i, j])^2;
      sO2 := sO2 + observed[i, j]^2;
      sOO := sOO + (observed[i, j] - meanO)^2;
      sMM := sMM + (simulated[i, j] - meanM)^2;
      sOM := sOM + (observed[i, j] - meanO) * (simulated[i, j] - meanM);
    end for;
  end for;
  chi2 := sumSquares / max(noiseSigma, 1.0e-300)^2;
  explained := if sO2 > 0.0 then 1.0 - sumSquares / sO2 else 0.0;
  correlation := if sOO > 0.0 and sMM > 0.0 then sOM / sqrt(sOO * sMM) else 0.0;
end imageResidual;

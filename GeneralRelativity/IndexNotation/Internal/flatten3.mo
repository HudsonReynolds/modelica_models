within GeneralRelativity.IndexNotation.Internal;
function flatten3 "Row-major flat components of a 4x4x4 array"
  input Real A[4, 4, 4];
  output Real flat[64];
algorithm
  for i in 1:4 loop
    for j in 1:4 loop
      for k in 1:4 loop
        flat[((i - 1) * 4 + (j - 1)) * 4 + k] := A[i, j, k];
      end for;
    end for;
  end for;
end flatten3;

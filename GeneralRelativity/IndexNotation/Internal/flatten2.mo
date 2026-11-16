within GeneralRelativity.IndexNotation.Internal;
function flatten2 "Row-major flat components of a 4x4 array"
  input Real A[4, 4];
  output Real flat[16];
algorithm
  for i in 1:4 loop
    for j in 1:4 loop
      flat[(i - 1) * 4 + j] := A[i, j];
    end for;
  end for;
end flatten2;

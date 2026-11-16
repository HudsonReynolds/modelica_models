within GeneralRelativity.IndexNotation.Internal;
function flatten4 "Row-major flat components of a 4x4x4x4 array"
  input Real A[4, 4, 4, 4];
  output Real flat[256];
algorithm
  for i in 1:4 loop
    for j in 1:4 loop
      for k in 1:4 loop
        for l in 1:4 loop
          flat[(((i - 1) * 4 + (j - 1)) * 4 + (k - 1)) * 4 + l] := A[i, j, k, l];
        end for;
      end for;
    end for;
  end for;
end flatten4;

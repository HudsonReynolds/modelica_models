within GeneralRelativity.Cosmology;
function tableValue "Linear interpolation in a table uniform in ln(eta), clamped at its ends"
  input Real lnEtaStart "ln(eta) of the first entry";
  input Real dLnEta "Spacing in ln(eta)";
  input Real values[:];
  input Real eta;
  output Real v;
protected
  Real u = (log(eta) - lnEtaStart) / dLnEta;
  Integer i;
algorithm
  if u <= 0.0 then
    v := values[1];
  elseif u >= size(values, 1) - 1 then
    v := values[size(values, 1)];
  else
    i := integer(u);
    v := values[i + 1] + (u - i) * (values[i + 2] - values[i + 1]);
  end if;
  annotation(Inline = true);
end tableValue;

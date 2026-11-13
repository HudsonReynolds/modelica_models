within GeneralRelativity.NewtonCartan;
function clockForm "Absolute-time one-form tau_a = d_a(ct) = {1, 0, 0, 0}"
  output Real tau[4] "Covariant components";
algorithm
  tau := {1.0, 0.0, 0.0, 0.0};
  annotation(Inline = true);
end clockForm;

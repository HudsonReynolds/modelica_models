within GeneralRelativity.Symmetry;
function timeTranslation "Time translation d_0 = d/d(ct); Killing for stationary metrics"
  extends GeneralRelativity.Interfaces.VectorFieldFunction;
algorithm
  xi := {1.0, 0.0, 0.0, 0.0};
  annotation(Documentation(info="<html><p>Parameters <code>p</code> are ignored.</p></html>"));
end timeTranslation;

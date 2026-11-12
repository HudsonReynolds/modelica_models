within GeneralRelativity;
package Constants "Physical constants and body parameters in SI units"
  constant Real c(unit = "m/s") = 299792458.0 "Speed of light";
  constant Real G(unit = "m3/(kg.s2)") = 6.67430e-11 "Newtonian gravitational constant";
  constant Real pi = 3.141592653589793;
  constant Real kappa(unit = "s2/(kg.m)") = 8.0 * pi * G / c^4 "Einstein coupling 8 pi G / c^4";

  constant Real GM_sun(unit = "m3/s2") = 1.32712440018e20 "Heliocentric gravitational constant";
  constant Real GM_earth(unit = "m3/s2") = 3.986004418e14 "Geocentric gravitational constant";
  constant Real GM_moon(unit = "m3/s2") = 4.9028e12 "Selenocentric gravitational constant";
  constant Real R_earth(unit = "m") = 6378137.0 "Earth equatorial radius";
  constant Real J2_earth = 1.08263e-3 "Earth second zonal harmonic";
  constant Real omega_earth(unit = "rad/s") = 7.2921150e-5 "Earth rotation rate";
  constant Real J_earth(unit = "kg.m2/s") = 5.86e33 "Earth spin angular momentum";
  constant Real earthMoonDistance(unit = "m") = 384400.0e3 "Mean Earth-Moon distance";
  constant Real omega_moon(unit = "rad/s") = 2.6617e-6 "Mean sidereal lunar orbital rate";

  function schwarzschildRadius "r_s = 2 G M / c^2"
    input Real GM(unit = "m3/s2") "Gravitational parameter G*M";
    input Real c(unit = "m/s") = Constants.c "Speed of light";
    output Real rs(unit = "m");
  algorithm
    rs := 2.0 * GM / c^2;
    annotation(Inline = true);
  end schwarzschildRadius;

  function geometricMass "Mass as a length, M = G M / c^2"
    input Real GM(unit = "m3/s2") "Gravitational parameter G*M";
    input Real c(unit = "m/s") = Constants.c "Speed of light";
    output Real M(unit = "m");
  algorithm
    M := GM / c^2;
    annotation(Inline = true);
  end geometricMass;
  annotation(Documentation(info="<html>
<p>SI values from CODATA 2018 and IERS Conventions 2010. The Sun's Schwarzschild
radius is <code>schwarzschildRadius(GM_sun)</code> = 2953 m; Earth's is 8.87 mm.
Geometric units set <code>c = G = 1</code> and measure mass in metres through
<code>geometricMass</code>.</p>
</html>"));
end Constants;

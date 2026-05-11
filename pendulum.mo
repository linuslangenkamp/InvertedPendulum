// Simple Inverted Pendulum with flat Modelica equations
// TODO: use standard components: see in Modelica.Mechanics.MultiBody.Examples.Elementary.Pendulum

model InvertedPendulum
  parameter Real M = 5;
  parameter Real m = 0.5;
  parameter Real l = 0.5;
  parameter Real g = 9.81;
  parameter Real d_cart = 0.15;
  parameter Real d_pend = 0.15;

  input Real tau;
  output Real s;
  output Real v;
  output Real phi;
  output Real vphi;

  Real a;
  Real alpha;

initial equation
  s = 0;
  v = 0;
  phi = 0.75 * Modelica.Constants.pi / 2;
  vphi = 0;

equation
  der(s) = v;
  der(phi) = vphi;

  a = (tau + m * sin(phi) * (l * vphi^2 + g * cos(phi)) - d_cart * v) / (M + m * sin(phi)^2);
  alpha = (-tau * cos(phi) - m * l * vphi^2 * sin(phi) * cos(phi)
           - (M + m) * g * sin(phi) - d_pend * vphi) / (l * (M + m * sin(phi)^2));

  der(v) = a;
  der(vphi) = alpha;

end InvertedPendulum;

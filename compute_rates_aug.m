function dX = compute_rates_aug(t, X)
  % Splits X into [rv; Φ], computes their time‐derivatives, and re‑stacks.

  global mu J2 Re

  % unpack
  rv  = X(1:6);                 % [r; v]
  Phi = reshape(X(7:42),6,6);   % state‑transition matrix

  % 1) Nonlinear two‑body + J2 rates:
  r = rv(1:3);  v = rv(4:6);
  drv = compute_rates_rv_perturbed(t, rv);

  % 2) Build Jacobian F(r):
  F = compute_F_j2(r);

  % 3) Variational equation: dΦ/dt = F * Φ
  dPhi = F * Phi;

  % pack up
  dX = [drv; dPhi(:)];
end
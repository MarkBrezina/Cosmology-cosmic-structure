src = open('/home/claude/lab/followups.py').read()
exec(src.split("print(f\"\\nbest mimic")[0].replace('print(f"uncorrected', 'pass  # print(f"uncorrected'))
class Smooth:
    """Same mean history as the tuned law, but dark energy kept smooth in space (keys on the average density only)."""
    def __init__(self, lm): self.lm, self.f_grid = lm, lm.f_grid
    def f(self, z): return self.lm.f(z)
    def g(self, z): return np.zeros_like(np.asarray(z, float))
sm = Smooth(lm_m)
r = optimize.minimize_scalar(lambda Om: sum(chi2_parts(sm, Om).values()), bounds=(0.1, 0.6), method="bounded")
f, D = growth(sm, float(r.x))
print(f"average-density version: SN+BAO+flows Δχ² vs ΛCDM {r.fun - c_lcdm['SN+BAO+flows']:+.2f}, predicted flow slope {f*D/REF_FD:.2f}")

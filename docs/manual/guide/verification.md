# Verifying a calculation

Two different questions are worth keeping apart.

*Verification* asks whether the code does what it was designed to do.
*Validation* asks whether that design describes nature, and needs a reference
outside the code.

## What the package checks for you

The test suite covers XML round-tripping as a fixed point, ReacLib rates on a
$T_9$ grid, screening factors, weak-rate interpolation, right-hand-side values
and their baryon residual, full trajectories, the analytic Jacobian against
finite differences, the treatment of photons and leptons, solver dispatch, the
spin degeneracy, exclusion of species without nuclear data, and detection of a
positivity projection that manufactures mass.

```bash
pytest -q
```

## What you should check on your own calculation

**Conservation.** For a reaction set that conserves baryon number, the
instantaneous residual

$$\delta_A(t)=\sum_i A_i\,\frac{dY_i}{dt}$$

should stay at round-off, and $\sum_i A_iY_i$ should stay at its initial value.
Track them along the whole trajectory, not only at the start.

**Solver agreement.** Repeat the calculation with `bdf`, `radau` and `lsoda`.
They should agree far more closely than any physical effect you care about.

**Tolerance convergence.** Tighten `rtol` and confirm the endpoint stops
moving.

**Equilibrium.** Where the physics allows it, run the calculation to its
stationary state and compare with `solve_nse`. This is the strongest check
available, because the two share no machinery. Read {doc}`pitfalls` on reverse
rates first, or the comparison will measure the rate library rather than the
code.

## Golden files and what they do not prove

`tests/golden/` holds frozen outputs, and
`nucnetpy.validation.regression_summary` compares against them. These detect
accidental changes in behaviour. They are snapshots of this implementation, so
they are not independent evidence that an answer is right, and they must not be
described as agreement with the original C++ code.

`validation/generate_golden.py` writes the reference records in the format the
suite consumes, so outputs from a NucNet Tools build can be substituted for the
self-generated ones without changing any test code.

## Reproducing the published results

```bash
python validation/demonstration.py \
    --nuclides nuclides.xml --reactions reaction_data.xml --zone zone.xml \
    --outdir results/

python validation/benchmark_performance.py \
    --nuclides nuclides.xml --reactions reaction_data.xml \
    --json results/benchmark.json --latex results/table_benchmark.tex
```

Both record the interpreter, NumPy, SciPy and platform versions alongside their
output.

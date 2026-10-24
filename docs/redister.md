# REDISTER

REDISTER is by Randy Gladstone; the original source is dated June 27, 1991.
Sources are under `src/glow_los/rt/redister_src`. Preserve the original header.
The wrapper's BSD-2-Clause license does not replace upstream source terms.

## Underflow message

The repeated message is:

```text
Note: The following floating-point exceptions are signalling: IEEE_UNDERFLOW_FLAG
```

A diagnostic build with `-ffpe-trap=underflow -g -fbacktrace` locates underflow
in REDISTER's `COEFF` calculation for the daytime example. Very small transfer
coefficients reach the floating-point underflow range. This is distinct from
an invalid operation, division by zero, or overflow.

The builder uses `-ffpe-summary=invalid,zero,overflow` for REDISTER. This changes
the exit summary only; it does not alter the calculation or enable fast math.
Serious floating-point diagnostics remain enabled. Rebuild with:

```bash
glow-los-build
```

For debugging, restore all summaries:

```bash
FFLAGS='-ffpe-summary=all' glow-los-build
```

See [GNU Fortran diagnostic options](https://gcc.gnu.org/onlinedocs/gfortran/Debugging-Options.html).
This addresses the reproduced message; unexpected nonfinite radiances still
need investigation for the specific atmosphere and geometry.

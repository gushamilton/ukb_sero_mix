#!/usr/bin/env python3
# Runtime monkey‑patch for pandas ≥ 2.0
import sys, pandas as pd
_orig_read_csv = pd.read_csv
def _patched_read_csv(filepath_or_buffer, *args, **kw):
    if args:                   # legacy positional separator
        kw.setdefault("sep", args[0])
        args = args[1:]
    return _orig_read_csv(filepath_or_buffer, *args, **kw)
pd.read_csv = _patched_read_csv

# hand execution over to the real quickdraws entry‑point
from quickdraws.scripts.quickdraws_step2 import main
sys.exit(main())


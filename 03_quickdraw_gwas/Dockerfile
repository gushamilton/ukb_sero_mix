# ------------------------------------------------------------------
#  Quickdraws Dockerfile — fixes legacy positional '\t' read_csv
#  and fails the build if *any* instance remains.
# ------------------------------------------------------------------
    FROM nvidia/cuda:11.8.0-runtime-ubuntu22.04

    # ---------- 1.  System packages ----------
    RUN apt-get update && apt-get install -y --no-install-recommends \
            build-essential \
            python3 \
            python3-pip \
            python3-dev \
            python-is-python3 \
        && rm -rf /var/lib/apt/lists/*
    
    # ---------- 2.  Latest pip tooling ----------
    RUN pip install --upgrade pip setuptools wheel
    
    # ---------- 3.  Install quickdraws + modern numpy/pandas ----------
    RUN pip install --no-cache-dir numpy pandas quickdraws
    
    # ---------- 4.  Auto‑patch every offending read_csv call ----------
    RUN python3 - <<'PY'
    import site, pathlib, re, sys, textwrap
    
    # Regex matches: pd.read_csv(<anything>, '\t')
    PAT = re.compile(r"pd\.read_csv\(\s*([^,]+?)\s*,\s*'\\\\t'\s*\)")
    
    patched = 0
    for sp_dir in site.getsitepackages():
        qroot = pathlib.Path(sp_dir, "quickdraws")
        if not qroot.exists():
            continue
        for py in qroot.rglob("*.py"):
            code = py.read_text()
            if PAT.search(code):
                new_code = PAT.sub(r"pd.read_csv(\1, sep='\\t')", code)
                py.write_text(new_code)
                patched += 1
    
    print(f"🔧  Patched {patched} positional read_csv call(s).")
    
    # Abort if *any* instance still lurks
    leftover = []
    for sp_dir in site.getsitepackages():
        qroot = pathlib.Path(sp_dir, "quickdraws")
        if not qroot.exists():
            continue
        for py in qroot.rglob("*.py"):
            if PAT.search(py.read_text()):
                leftover.append(str(py))
    if leftover:
        sys.stderr.write(textwrap.dedent(f"""
            ❌  Unpatched positional read_csv calls remain:
            {chr(10).join(leftover)}
            """))
        sys.exit(1)
    
    print("✅  All occurrences patched.")
    PY
    
    # ---------- 5.  Import‑time sanity check ----------
    RUN python3 - <<'PY'
    import pandas, numpy, importlib
    print("pandas:", pandas.__version__, "numpy:", numpy.__version__)
    importlib.import_module("quickdraws.scripts.quickdraws_step2")
    print("✅  quickdraws.step2 imports cleanly with patched read_csv.")
    PY
    
    # ---------- 6.  Image metadata ----------
    WORKDIR /app
    CMD ["/bin/bash"]
    
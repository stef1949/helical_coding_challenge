# Helical2 Project Notes

## Installing the Helical SDK on Windows

The upstream `helical` PyPI wheel depends on `pytorch==2.6.0 --index-url https://download.pytorch.org/whl/cu126` & `louvain==0.8.2`.  Building from source will fail on Windows because it requires
native igraph/CMake tooling that is not readily available.

Use the helper script instead:

```powershell
.\.venv\Scripts\python.exe scripts/install_helical.py
```

The script patches the official Helical wheel to drop the `louvain` requirement
and then installs it together with the rest of Helical's dependencies.  Re-run
it any time you recreate the virtual environment.

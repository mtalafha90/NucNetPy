"""Allow the command line to be invoked as ``python -m nucnetpy.cli``.

The console script installed as ``nucnetpy`` is the usual entry point, but it
requires the environment's bin directory to be on PATH. Running the module
works from any interpreter that can import the package, which is what tests and
scripts want.
"""
from . import main

raise SystemExit(main())

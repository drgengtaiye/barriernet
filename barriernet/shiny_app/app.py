"""barriernet web app (Shiny for Python).

Run it here with `shiny run app.py`, or let Shiny Server 1.5.21 or later serve a folder made by make_bundle.sh, which
puts a copy of the barriernet package next to this file. Either way the package is found without installing it."""
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
for root in (os.path.dirname(HERE), HERE):              # the repository layout, then a bundle made by make_bundle.sh
    if os.path.isfile(os.path.join(root, "barriernet", "__init__.py")) and root not in sys.path:
        sys.path.insert(0, root)
os.environ.setdefault("MPLCONFIGDIR", os.path.join(tempfile.gettempdir(), "barriernet-matplotlib"))

from barriernet.app import app  # noqa: E402,F401

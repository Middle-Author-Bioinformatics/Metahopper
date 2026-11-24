"""Shared report browsing controls; also embedded by the signed-in viewer."""
from pathlib import Path

BASE = Path(__file__).parent
CONTROL_CSS = (BASE / "metahopper_report_controls.css").read_text()
CONTROL_JS = (BASE / "metahopper_report_controls.js").read_text()

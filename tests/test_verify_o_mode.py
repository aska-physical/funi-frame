"""`python3 -O` strips every assert statement. These tests pin the verifier's checks against that.

Run from the repository root:  python3 -B -m unittest discover -s tests
"""
import ast
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# Make must-fail control 1 fail, then build the certificate.
BREAK_ONE_CONTROL = f"""
import sys
sys.path.insert(0, {str(ROOT)!r})
import verify
roof = verify.solve_funicular_roof()
verify.solve_funicular_roof = lambda: {{**roof, "span_over_delta_vault": 100.0}}
try:
    verify.build_certificate()
except RuntimeError as error:
    print("STOPPED", error)
else:
    print("NOT STOPPED")
"""


def run(*args):
    return subprocess.run([sys.executable, "-B", *args], cwd=ROOT, capture_output=True, text=True)


class OptimiseMode(unittest.TestCase):
    def test_verifier_has_no_assert_statements(self):
        tree = ast.parse((ROOT / "verify.py").read_text(encoding="utf-8"))
        self.assertEqual([n.lineno for n in ast.walk(tree) if isinstance(n, ast.Assert)], [])

    def test_clean_tree_passes_with_and_without_optimise(self):
        for flags in ([], ["-O"]):
            with self.subTest(flags=flags):
                result = run(*flags, "verify.py")
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn("[PASS]", result.stdout)

    def test_failing_control_stops_the_run_with_and_without_optimise(self):
        for flags in ([], ["-O"]):
            with self.subTest(flags=flags):
                result = run(*flags, "-c", BREAK_ONE_CONTROL)
                self.assertIn("STOPPED", result.stdout, result.stderr)
                self.assertIn("must_fail_01", result.stdout)


if __name__ == "__main__":
    unittest.main()

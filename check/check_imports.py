import platform
import sys
import subprocess

required_packages = {
    "MDAnalysis": "MDAnalysis",
    "psfgen": "psfgen",
    "numpy": "numpy",
    "vmd": "vmd",
    "PyQt5": "PyQt5"
}

missing = []

# Check Python packages
for pkg, import_name in required_packages.items():
    try:
        __import__(import_name)
        print(f"✔️  {pkg} imported successfully.")
    except ImportError:
        print(f"❌ {pkg} NOT installed or import failed.")
        missing.append(pkg)

# Check AmberTools (tleap)
def check_tleap():
    try:
        result = subprocess.run(['tleap', '-h'], capture_output=True, text=True)
        if result.returncode == 0 or 'usage:' in result.stdout.lower():
            print("✔️  AmberTools is installed and accessible.")
            return True
        else:
            print("❌ AmberTools might not be properly installed.")
            return False
    except FileNotFoundError:
        print("❌ AmberTools is NOT installed or not in PATH.")
        return False

ambertools_ok = check_tleap()
if not ambertools_ok:
    missing.append("AmberTools")


if not missing:
    print("ALL_OK")
    sys.exit(0)
else:
    print(" ".join(missing))  # <-- Output just package names
    sys.exit(1)

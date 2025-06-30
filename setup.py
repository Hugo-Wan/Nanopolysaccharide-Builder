from setuptools import setup, find_packages
import subprocess
import json
import os
import atexit
import sys
import platform
import shutil
import requests
import zipfile
import shutil
import stat

def find_conda():
    """ Attempt to find a 'conda' executable in the system path. """
    conda_exec = 'conda'  # Default, works if 'conda' is in the PATH
    for path in os.environ.get('PATH', '').split(os.pathsep):
        potential_conda = os.path.join(path, conda_exec)
        if os.path.exists(potential_conda) and os.access(potential_conda, os.X_OK):
            return potential_conda
    return None  # Return None if not found


def make_all_files_executable(path):
    for root, dirs, files in os.walk(path):
        for f in files:
            file_path = os.path.join(root, f)
            try:
                # Add user read/write/execute permission
                os.chmod(file_path, os.stat(file_path).st_mode | stat.S_IRWXU)
            except Exception as e:
                print(f"Warning: Could not change permission for {file_path}: {e}")


def get_missing_packages():
    # Check current directory and one level up for check_imports.py
    possible_paths = [
        "./check/check_imports.py"
    ]
    script_path = None
    for path in possible_paths:
        if os.path.isfile(path):
            script_path = path
            break

    if script_path is None:
        print("check_imports.py not found in current or parent directory!")
        sys.exit(1)

    # Run check_imports.py and get output
    result = subprocess.run([sys.executable, script_path], capture_output=True, text=True)
    last_line = result.stdout.strip().split("\n")[-1]
    if last_line == "ALL_OK":

        return []
        
    else:
        conda_path = find_conda()
        print("Cleaning conda packages and tarballs, please wait...")
        subprocess.run([conda_path, "clean", "--packages", "--tarballs", "--yes"], check=True, capture_output=True)
        return last_line.split()


def attempt_install():
    sys_type = platform.system()
    machine = platform.machine()
    conda_path = find_conda()
    print(f"System type : {sys_type}")
    if sys_type == "Linux":
        print(f"Installing NPB in {machine} system")
        
        if not conda_path:
            print("Conda executable not found. Please ensure Conda is installed and in the system PATH.")
            return
        missing = get_missing_packages()
        if missing:
            print("Missing packages to install:", missing)
            try:
                if "AmberTools" in missing:
                    print(f"Installing ambertools, please wait...")
                    subprocess.run([conda_path, "install", "conda-forge::ambertools", "-y"], check=True, capture_output=True)
                if "vmd" in missing or "vmd-python" in missing:
                    print(f"Installing vmd-python, please wait...")
                    subprocess.run([conda_path, "install", "-c", "conda-forge", "vmd-python", "-y"], check=True, capture_output=True)
                if "vmd" in missing or "vmd-python" in missing:    
                    print(f"Installing psfgen-python, please wait...")
                    subprocess.run([conda_path, "install", "-c", "conda-forge", "psfgen", "-y"], check=True, capture_output=True)
                if "PyQt5" in missing:     
                    print("Installing PyQt5, please wait...")
                    subprocess.run([conda_path, "install", "-c", "conda-forge", "pyqt", "-y"], check=True, capture_output=True)
                if "MDAnalysis" in missing:    
                    print(f"Installing mdanalysis, please wait...")
                    subprocess.run([conda_path, "install", "-c", "conda-forge", "mdanalysis", "-y"], check=True, capture_output=True)
                if "numpy" in missing:    
                    print(f"Installing numpy, please wait...")
                    subprocess.run([conda_path, "install", "-c", "conda-forge", "numpy", "-y"], check=True, capture_output=True)
                missing = get_missing_packages()  
                if missing:
                    print("Missing packages for automatic installation:", missing, "Please install them manually.")
                else:
                    print("\nRequired dependencies installed successfully.")
                    print("\n\nEnjoy using NPB !!!")     
                    
            except subprocess.CalledProcessError as e:
                print(f"\nFailed to install packages via Conda:\nSTDOUT: {e.stdout.decode()}\nSTDERR: {e.stderr.decode()}")

        else:
            print("\nAll packages are present!")
            print("\nEnjoy using NPB !!!")


            
    elif sys_type == "Darwin" and machine in ("arm64", "aarch64"):
        print(f"Installing NPB in Darwin (mac-osx-arm) system")
        
        if not conda_path:
            print("Conda executable not found. Please ensure Conda is installed and in the system PATH.")
            return
        
        missing = get_missing_packages()
        if missing:
            print("Missing packages to install:", missing)
            try:
                if "AmberTools" in missing:
                    print(f"Installing ambertools, please wait...")
                    subprocess.run([conda_path, "install", "conda-forge::ambertools", "-y"], check=False, capture_output=True)    
                if "MDAnalysis" in missing:
                    print("Installing mdanalysis, please wait...")
                    subprocess.run([conda_path, "run", "pip", "install", "--upgrade", "MDAnalysis"], check=True, capture_output=True)
                if "numpy" in missing:
                    print(f"Installing numpy, please wait...")
                    subprocess.run([conda_path, "install", "-c", "conda-forge", "numpy", "-y"], check=True, capture_output=True)
                if "PyQt5" in missing:     
                    print("Installing PyQt5, please wait...")
                    subprocess.run([conda_path, "run", "pip", "install", "PyQt5"], check=True, capture_output=True)

            
                # --- psfgen-python ---
                if "psfgen" in missing:
                    psfgen_url = "https://github.com/Eigenstate/psfgen/archive/refs/heads/master.zip"
                    psfgen_zip_path = "psfgen-python.zip"
                    psfgen_extract_dir = "psfgen-master"
            
                    print("Downloading psfgen-python from GitHub...")
                    with requests.get(psfgen_url, stream=True) as r:
                        with open(psfgen_zip_path, 'wb') as f:
                            shutil.copyfileobj(r.raw, f)
                    print("Download complete.")
            
                    print("Unzipping psfgen-python...")
                    with zipfile.ZipFile(psfgen_zip_path, 'r') as zip_ref:
                        zip_ref.extractall(".")
                    make_all_files_executable(psfgen_extract_dir)
            
                    if os.path.isdir(psfgen_extract_dir):
                        os.chdir(psfgen_extract_dir)
                        print("Building and installing psfgen-python...")
                        try:
                            subprocess.run([sys.executable, "setup.py", "install"], check=True)
                            print("psfgen-python installed successfully.")
                            
                        except subprocess.CalledProcessError:
                            print("Failed to build/install psfgen-python.")
                        os.chdir("..")
                    else:
                        print(f"Could not find extracted directory: {psfgen_extract_dir}")
                    os.remove(psfgen_zip_path)
                

                if "vmd" in missing:
                    print("Installing PyQt5, please wait...")
                    subprocess.run([conda_path, "install", "conda-forge::vmd-python", "-y"], check=True, capture_output=True) 
                ##check if conda install not successfully
                missing = get_missing_packages() 
                if "vmd" in missing:  
                    conda_path = shutil.which("conda")
                    if conda_path is None:
                        raise RuntimeError("Conda not found in PATH.")
                    
                    print("Installing expat, please wait...")
                    subprocess.run([conda_path, "install", "-c", "conda-forge", "expat=2.7.0", "-y"], check=True, capture_output=True)
                    print("Installing libnetcdf, please wait...")
                    subprocess.run([conda_path, "install", "conda-forge::libnetcdf", "-y"], check=True, capture_output=True)
                    
                    conda_lib_path = os.path.join(os.environ.get('CONDA_PREFIX', ''), 'lib')
                    bashrc_path = os.path.expanduser("~/.bash_profile")
                    export_line = f'\n# Added by setup script for expat/libnetcdf\nexport DYLD_LIBRARY_PATH="{conda_lib_path}:$DYLD_LIBRARY_PATH"\n'
                    with open(bashrc_path, 'a') as bashrc:
                        bashrc.write(export_line)
                    
                    if shutil.which("swig") is None:
                        print("SWIG not found. Attempting to install with Homebrew...")
                        try:
                            subprocess.run(["brew", "install", "swig"], check=True)
                            print("SWIG installed successfully via Homebrew.")
                        except subprocess.CalledProcessError:
                            print("Failed to install SWIG. Please install it manually.")
                    
                    # --- vmd-python ---
                    vmd_python_url = "https://github.com/Eigenstate/vmd-python/archive/refs/heads/master.zip"
                    vmd_zip_path = "vmd-python.zip"
                    vmd_extract_dir = "vmd-python-master"
                    
                    print("Downloading vmd-python from GitHub...")
                    with requests.get(vmd_python_url, stream=True) as r:
                        with open(vmd_zip_path, 'wb') as f:
                            shutil.copyfileobj(r.raw, f)
                    print("Download complete.")
                    
                    print("Unzipping vmd-python...")
                    with zipfile.ZipFile(vmd_zip_path, 'r') as zip_ref:
                        zip_ref.extractall(".")
                    make_all_files_executable(vmd_extract_dir)
                    
                    if os.path.isdir(vmd_extract_dir):
                        os.chdir(vmd_extract_dir)
                        print("Building vmd-python...")
                        try:
                            subprocess.run([sys.executable, "setup.py", "build"], check=True)
                            subprocess.run([sys.executable, "setup.py", "install"], check=True)
                            print("vmd-python installed successfully.")
                        except subprocess.CalledProcessError:
                            print("Failed to build/install vmd-python.")
                        os.chdir("..")
                    else:
                        print(f"Could not find extracted directory: {vmd_extract_dir}")
                    os.remove(vmd_zip_path)
                missing = get_missing_packages()  
                if missing:
                    print("Missing packages for automatic installation:", missing, "Please install them manually.")
                else:
                    print("\nRequired dependencies installed successfully.")
                    print("\n\nEnjoy using NPB !!!")   

            except subprocess.CalledProcessError as e:
                print(f"\nFailed to install all the required dependencies:\nSTDOUT: {e.stdout.decode()}\nSTDERR: {e.stderr.decode()}")

        else:
            print("\nAll packages are present!")
            print("\nEnjoy using NPB !!!")

    elif sys_type == "Darwin" and machine == "x86_64":
        print(f"Installing NPB in Darwin (mac-osx) system")
        if not conda_path:
            print("Conda executable not found. Please ensure Conda is installed and in the system PATH.")
            return
        
        missing = get_missing_packages()
        if missing:
            print("Missing packages to install:", missing)
            try:
    
                if "AmberTools" in missing:
                    print(f"Installing ambertools, please wait...")
                    subprocess.run([conda_path, "install", "conda-forge::ambertools", "-y"], check=True, capture_output=True)    
                if "MDAnalysis" in missing:
                    print("Installing mdanalysis, please wait...")
                    subprocess.run([conda_path, "run", "pip", "install", "--upgrade", "MDAnalysis"], check=True, capture_output=True)
                if "numpy" in missing:
                    print(f"Installing numpy, please wait...")
                    subprocess.run([conda_path, "install", "-c", "conda-forge", "numpy", "-y"], check=True, capture_output=True)
                if "PyQt5" in missing:    
                    print("Installing PyQt5, please wait...")
                    subprocess.run([conda_path, "run", "pip", "install", "PyQt5"], check=True, capture_output=True)
                if "vmd" in missing or "vmd-python" in missing:
                    print(f"Installing vmd-python, please wait...")
                    subprocess.run([conda_path, "install", "-c", "conda-forge", "vmd-python", "-y"], check=True, capture_output=True)
                if "psfgen" in missing:
                    print(f"Installing vmd-python, please wait...")
                    subprocess.run([conda_path, "install", "-c", "conda-forge::psfgen", "-y"], check=True, capture_output=True)   
                #check if not install properly for psfgen
                missing = get_missing_packages()                 
                if "psfgen" in missing:
                    psfgen_url = "https://github.com/Eigenstate/psfgen/archive/refs/heads/master.zip"
                    psfgen_zip_path = "psfgen-python.zip"
                    psfgen_extract_dir = "psfgen-master"
            
                    print("Downloading psfgen-python from GitHub...")
                    with requests.get(psfgen_url, stream=True) as r:
                        with open(psfgen_zip_path, 'wb') as f:
                            shutil.copyfileobj(r.raw, f)
                    print("Download complete.")
            
                    print("Unzipping psfgen-python...")
                    with zipfile.ZipFile(psfgen_zip_path, 'r') as zip_ref:
                        zip_ref.extractall(".")
                    make_all_files_executable(psfgen_extract_dir)
            
                    if os.path.isdir(psfgen_extract_dir):
                        os.chdir(psfgen_extract_dir)
                        print("Building and installing psfgen-python...")
                        try:
                            subprocess.run([sys.executable, "setup.py", "install"], check=True)
                            print("psfgen-python installed successfully.")
                            
                        except subprocess.CalledProcessError:
                            print("Failed to build/install psfgen-python.")
                        os.chdir("..")
                    else:
                        print(f"Could not find extracted directory: {psfgen_extract_dir}")
                    os.remove(psfgen_zip_path)
                missing = get_missing_packages()  
                if missing:
                    print("Missing packages for automatic installation:", missing, "Please install them manually.")
                else:
                    print("\nRequired dependencies installed successfully.")
                    print("\n\nEnjoy using NPB !!!")     
            except subprocess.CalledProcessError as e:
                print(f"\nFailed to install all the required dependencies:\nSTDOUT: {e.stdout.decode()}\nSTDERR: {e.stderr.decode()}")
        else:
            print("\nAll packages are present!")
            print("\nEnjoy using NPB !!!")


def get_install_dir():
    """Get the installation directory safely."""
    if '__file__' in globals():
        return os.path.dirname(os.path.abspath(__file__))
    else:
        return os.getcwd()

def create_config():
    """Create or update the config.json file based on installation directory."""
    install_dir = get_install_dir()
    config = {
        "main_folder_path": install_dir,
        "enable_logging": True,
        "log_level": "INFO"
    }
    config_path = os.path.join(install_dir, 'config.json')
    with open(config_path, 'w') as f:
        json.dump(config, f, indent=4)
    print("\nConfiguration file created at:", config_path)
    
if __name__ == "__main__":
    answer = input(
    "Would you like to attempt automatic required package installation? [y/N]:\n"
    "-------------------------------------------------------------\n"
    "|If you have installed all dependencies, you can just type N|\n"
    "-------------------------------------------------------------\n"
    ).strip().lower()

    if answer == "y":
        attempt_install()
    else:
        print("Skipping required dependencies install step.")
    atexit.register(create_config)  # Register the config file creation to run on exit.

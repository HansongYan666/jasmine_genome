import shutil
import sys
import os

def get_tool_paths(software_name):
    software_name_test = shutil.which(software_name)
    if not software_name_test:
        sys.exit('Error: could not find '+software_name)

def get_jar_path(jar_path):
    if jar_path and os.path.isfile(jar_path):
        return os.path.abspath(jar_path)
    for directory in os.environ['PATH'].split(':'):
        try:
            path_files = [f for f in os.listdir(directory)
                          if os.path.isfile(os.path.join(directory, f))]
        except FileNotFoundError:
            path_files = []
        pilon_jars = [f for f in path_files if f.startswith(jar_path) and f.endswith('.jar')]
        if pilon_jars:
            return os.path.join(directory, sorted(pilon_jars)[-1])  # return the latest version
    return None
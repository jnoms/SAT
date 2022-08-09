import pathlib
import os
import sys


# ------------------------------------------------------------------------------------ #
# Misc
# ------------------------------------------------------------------------------------ #
def talk_to_me(msg):
    """
    Prints message to the screen, but prefixes it with the parent pyscript name.
    """
    base = os.path.basename(sys.argv[0])
    print(f"{base}: {msg}")


# ------------------------------------------------------------------------------------ #
# Read/write
# ------------------------------------------------------------------------------------ #
def make_output_dir(path, is_dir=False):
    """
    Finds the dirname of the path, and makes the directory if it doesn't exist.
    Notably, if the path is the a directory it skips the dirname finding and just makes
    that directory
    """
    if is_dir == False:
        out_dir = os.path.dirname(path)
        pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)
    else:
        pathlib.Path(path).mkdir(parents=True, exist_ok=True)


if __name__ == "__main__":
    msg = "This script has utilities and functions. Don't call it directly!"
    raise ValueError(msg)

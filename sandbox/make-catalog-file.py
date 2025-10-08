import glob
import sys

location = sys.argv[1]

files = glob.glob(f"{location}/*.hdf5")

if len(files) > 0:

    for file in files:
        name = file.split("-")[3]
        print(f"- name: {name}")
        print(f"  metafile: {file}")
        print(f"  analysis: C00:Mixed")

else:
        
    files = glob.glob(f"{location}/*_cosmo.h5")
    for file in files:
        name = "_".join(file.split("-")[3].split("_")[:2])
        print(f"- name: {name}")
        print(f"  metafile: {file}")
        print(f"  analysis: C01:Mixed")

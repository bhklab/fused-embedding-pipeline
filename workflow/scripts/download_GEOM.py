import sys
import pandas as pd
import requests
import urllib
from damply import dirs
import yaml
import os 
import pathlib
from pathlib import Path


# Open the YAML file in file sources
with open(dirs.CONFIG / "downloads.yaml", 'r') as file:
    # Load the contents safely
    data = yaml.safe_load(file)


data = data['GEOM']
base_dir= dirs.RAWDATA / data['dir']
os.makedirs(base_dir, exist_ok=True)
data.pop("dir")

for k in data.keys():
    link = data[k]
    file_name = link.split("/")[-1].split("?")[0]
    urllib.request.urlretrieve(link, base_dir/file_name)

# config

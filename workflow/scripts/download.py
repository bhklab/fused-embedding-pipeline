import urllib.request
from pathlib import Path

import yaml
from damply import dirs

with open(dirs.CONFIG / "downloads.yaml") as file:
    config = yaml.safe_load(file)

inputs = {
    "HDD": ["colData"],
    "JUMP": ["colData", "cpcnn_parquet"],
    "GEOM": ["hp_parquet"],
    "LINCS": ["signatures"],
}

for source, keys in inputs.items():
    directory = dirs.RAWDATA / config[source]["dir"]
    directory.mkdir(parents=True, exist_ok=True)
    for key in keys:
        url = config[source][key]
        destination = directory / url.split("/")[-1].split("?")[0]
        if destination.exists():
            print(f"Already downloaded: {destination}")
            continue
        print(f"Downloading {source}: {key}", flush=True)
        partial = Path(str(destination) + ".part")
        urllib.request.urlretrieve(url, partial)
        partial.replace(destination)

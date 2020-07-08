# import shutil
# import urllib.request as request
# from contextlib import closing

# URL = "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/knownCanonical.txt.gz"

# with closing(request.urlopen(URL)) as uh:
#     with open("knownCanonical.txt.gz", "wb") as fh:
#         shutil.copyfileobj(uh, fh)


import requests, sys

server = "https://rest.ensembl.org"
ext = "/vep/human/id/COSM476?"

r = requests.get(server + ext, headers={"Content-Type": "application/json"})

if not r.ok:
    r.raise_for_status()
    sys.exit()

decoded = r.json()
print(repr(decoded))


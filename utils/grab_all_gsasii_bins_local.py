#!/usr/bin/env /Users/y8z/opt/anaconda3/envs/g2python311/bin/python

import urllib.request
import tarfile
import os

thetarfile = "https://github.com/GSASII/binarytest/releases/download/v1.0.1/mac_64_p3.11_n1.26.tgz"

os.chdir("/Users/y8z/Dev/gsasii_zyp/bindist")

ftpstream = urllib.request.urlopen(thetarfile)
thetarfile = tarfile.open(fileobj=ftpstream, mode="r|gz")
thetarfile.extractall()

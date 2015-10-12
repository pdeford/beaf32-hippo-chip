import subprocess

f = open("~/beaf32-hippo-chip/beaf_links.txt")
for line in f:
    subprocess.call("curl -O %s" % line, shell=True)
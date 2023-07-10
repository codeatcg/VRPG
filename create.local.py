
#!/usr/bin/env python
import re

headFile = "templates/vrpg/local.head"
indexFile = "templates/vrpg/index.html"
localFile = "templates/vrpg/local.index.html"

headStr = ""
pat1 = re.compile('\s*<head>')
pat2 = re.compile('\s*<\/head>')
with open(headFile) as fh:
    headStr = fh.read()
    
with open(indexFile) as dxh,open(localFile,'w') as wfh:
    flag1 = False
    flag2 = False
    for line in dxh:
        if not flag1:
            wfh.write(line)
            if pat1.search(line):
                wfh.write(headStr)
                flag1 = True
        else:
            if flag2:
                wfh.write(line)
            else:
                if pat2.search(line):
                    wfh.write(line)
                    flag2 = True



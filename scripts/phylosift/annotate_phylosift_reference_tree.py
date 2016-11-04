#!/usr/bin/python
import os, sys
#refpkgdir = sys.argv[1]
#jplacedir = sys.argv[2]
refpkgdir = os.environ['refpkg']
jplacedir = os.environ['resjplace']
markertag = os.environ['marker']

# copy annotation from a jplace file (from phylosift pipeline)
lnfjplace = os.listdir(jplacedir)
nfjplace = None
for nf in lnfjplace:
  if nf.endswith('.jplace'):
    nfjplace = nf
    print "extract tree annotation from jplace file %s"%nfjplace
    break

finjp = open('%s/%s'%(jplacedir, nfjplace))
s = finjp.read()
finjp.close()
o = eval(s)
fouttree = open('%s/%s.annotated/annotatedreftreefromjplace.tree'%(refpkg, markertag), 'w')
fouttree.write( o['tree'])
fouttree.close()

# modify the CONTENTS.json index file
fincon = open('%s/%s.annotated/CONTENTS.json'%(refpkg, markertag), 'r')
lines = fincon.readlines()
fincon.close()
foutcon = open('%s/%s.annotated/CONTENTS.json'%(refpkg, markertag), 'w')
key = None
for line in lines:
  lsr = line.strip(' \t\n')
  if lsr.endswith('{'): key = lsr.split('"')[1]
  if lsr.endswith('}'): key = None
  if key=="files" and lsr.split('"')[1]=="tree":
    lsp = line.split('"')
    lsp[3] = "annotatedreftreefromjplace.tree"
    foutcon.write('"'.join(lsp))
  else:
    foutcon.write(line)

foutcon.close()

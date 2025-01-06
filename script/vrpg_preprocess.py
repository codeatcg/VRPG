

import argparse
import re
import os
import gzip
import sys
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),"module"))

import minipg

parser = argparse.ArgumentParser()

parser.add_argument('--sep',help="delimiter between sample and haplotype names, by default: '#'",default="#")

parser.add_argument('--rGFA',help="Input rGFA file",default=None)
parser.add_argument('--gafList',help="Input GAF file list",default=None)
parser.add_argument('--minMQ',help="Minimal mapping quality that are allowed when process GAF files, by default: 0",type=int,default=0)
parser.add_argument('--ncov',help="Disable calculating node coverage for each assembly",action="store_true",default=False)
parser.add_argument('--outDir',help="Output directory",default=None,required=True)

parser.add_argument('--minigraph',help="The absolute file path of minigraph executable file, this is used when rGFA file does not exist",default=None)
parser.add_argument('--asmList',help="Input assembly file list. When rGFA file does not exist this option can be taken to build the reference pangenome graph",default=None)
parser.add_argument('--thread',help="Number of threads, by default: 10",type=int,default=10)
parser.add_argument('--graphOpt',help="minigraph options (except '-cxggs -t') to generate graph, e.g.: '-q 5 -l 100000 -L 50'",default=None)

parser.add_argument('--index',help="Index the graph to speed up data visualization",action="store_true",default=False)
parser.add_argument('--range',help="Range size for indexing the graph, by default: 2000",type=int,default=2000)
parser.add_argument('--xDep',help="Search depth when creating graph indexes, by default: 100",type=int,default=100)

paras = parser.parse_args()


def checkOut(outDir):
    if os.path.exists(outDir):
        print("Error: output directory has already existed! {}".format(outDir))
        exit(1)
    try:
        os.mkdir(outDir)
    except:
        print("Error: output directory create failed! {}".format(outDir))
                
def changeHeader(sep,asmList,outDir):  
    newDir = os.path.join(outDir,"newAsmDir")
    os.mkdir(newDir)
    
    pat = re.compile("\s+")
    pat2 = re.compile("^\s*$")
    nameList = []
    dotList = []
    suffix = []
    with open(asmList) as al:
        for line in al:
            if line.startswith('#'):
                continue
            if pat2.match(line):
                continue
            line = line.strip()
            arr = pat.split(line)
            fName = arr[0].replace(sep,'.')
            
            nameList.append(arr[0])
            dotList.append(fName)
            flag = False
            
            if arr[1].endswith(".gz"):
                try:
                    fh = gzip.open(arr[1],'rb')
                    flag = True
                    suffix.append(".fa.gz")
                except:
                    print("Error: assembly file open failed! {}".format(arr[1]))
            else:
                suffix.append(".fa")
                try:
                    fh = open(arr[1])
                except:
                    print("Error: assembly file open failed! {}".format(arr[1]))
            
            if flag:
                try:
                    outFile = os.path.join(newDir,fName + ".fa.gz")
                    tfh = gzip.open(outFile,"wb")
                except:
                    print("Error: output file open failed! {}".format(outFile))
            else:
                try:
                    outFile = os.path.join(newDir,fName + ".fa")
                    tfh = open(outFile,"w")
                except:
                    print("Error: output file open failed! {}".format(outFile))
            
            if flag:
                for faLine in fh:
                    faLine = bytes.decode(faLine)
                    if faLine.startswith('>'):
                        tarr = pat.split(faLine)
                        preName = tarr[0].replace('>','>' + arr[0] + sep)
                        tfh.write(str.encode(preName + "\n"))
                    else:
                        tfh.write(str.encode(faLine))
            else:
                for faLine in fh:
                    if faLine.startswith('>'):
                        tarr = pat.split(faLine)
                        preName = tarr[0].replace('>','>' + arr[0] + sep)
                        tfh.write(preName + "\n")
                    else:
                        tfh.write(faLine)
                
            tfh.close()
            fh.close()
    return {"nameList":nameList,"dotList":dotList,"suffix":suffix}    

def rgraph(dotList,suffix,minigraph,graphOpt,thread,outDir):
    newAsmDir = os.path.join(outDir,"newAsmDir")
    newAsmList = [ os.path.join(newAsmDir, x + y) for x,y in zip(dotList,suffix)]
    allFa = " ".join(newAsmList)
    
    graphDir = os.path.join(outDir,"newGraph")
    os.mkdir(graphDir)
    gfaFile = os.path.join(graphDir,"out.gfa")
    
    command = ""
    if graphOpt is not None:
        command = minigraph + " -cxggs -t " + thread + " " + graphOpt + " " + allFa + " > " + gfaFile
    else:
        command = minigraph + " -cxggs -t " + thread + " " + allFa + " > " + gfaFile
    try:
        os.system(command)
    except:
        print("Error: Graph create failed. Please check the options!")
    
def mapAsm(nameList,dotList,suffix,minigraph,thread,outDir):
    newAsmDir = os.path.join(outDir,"newAsmDir")
    newAsmList = [ os.path.join(newAsmDir, x + y) for x,y in zip(dotList,suffix)]
    
    mapDir = os.path.join(outDir,"mapDir")
    os.mkdir(mapDir)
    mapListFile = os.path.join(mapDir,"gaf.list")
    
    graphDir = os.path.join(outDir,"newGraph")
    gfaFile = os.path.join(graphDir,"out.gfa")
    
    with open(mapListFile,'w') as fh:
        #for faFile,name in zip(newAsmList,dotList):
        endPoint = len(dotList)
        for i in range(endPoint):
            outGAF = os.path.join(mapDir,dotList[i] + ".gaf")
            command = minigraph + " -cxasm -t " + thread + " --vc " + gfaFile + " " + newAsmList[i] + " -o " + outGAF
            os.system(command)
            
            fh.write(nameList[i] + "\t" + outGAF + "\n")

def reduceGFA(sep,gfaFile,nodeFile,edgeFile,comChrFile,ndAsmLFile):
    flag = False
    if gfaFile.endswith(".gz"):
        fh = gzip.open(gfaFile,"rb")
        flag = True
    else:
        fh = open(gfaFile)
        
    comChrSet = set()
    asmSet = set()
    
    with open(nodeFile,'w') as nfh, open(edgeFile,'w') as efh, open(comChrFile,'w') as cfh, open(ndAsmLFile,'w') as afh:
        nfh.write("#Segment\tChr\tStart\tEnd\tLen\tRefOrNot\n")
        efh.write("#Source\tTarget\tOrigin1\tOrigin2\n")
        for line in fh:
            if flag:
                line = bytes.decode(line)
            line = line.strip()
            arr = line.split("\t")
            if arr[0] == 'S':
                segment = arr[1].replace("s","")
                slen = arr[3].split(":")[2]
                schr = arr[4].split(":")[2]
                
                tAsm = ""
                if sep not in schr: 
                    tAsm = "REF" + sep + "0"
                    schr = tAsm + sep + schr
                else:
                    asmArr = schr.split(sep)
                    tAsm = asmArr[0] + sep + asmArr[1]
                    
                start = arr[5].split(":")[2]
                off_start = str(int(start) + 1)
                end = str(int(start) + int(slen))
                refOrNot = arr[6].split(":")[2]
                nfh.write(segment + "\t" + schr + "\t" + off_start + "\t" + end +"\t" + slen + "\t" + refOrNot +"\n")
                
                if schr not in comChrSet:
                    comChrSet.add(schr)
                    
                if tAsm not in asmSet:
                    asmSet.add(tAsm)
                    afh.write(tAsm + "\n")
                    
            elif arr[0] == 'L':
                source = arr[1].replace("s","")
                origin1 = arr[2]
                target = arr[3].replace("s","")
                origin2 = arr[4]
                efh.write(source + "\t" + target + "\t" + origin1 + "\t" + origin2 + "\n")
                
        for tchr in sorted(comChrSet):
            cfh.write(tchr + "\n")
            
        #for asm in sorted(asmSet):
        #    afh.write(asm + "\n")
    
    fh.close()

  
def refChrList(sep,nodeFile,outChrList):
    with open(nodeFile) as nf,open(outChrList,'w') as cl:
        #cl.write("#Chr\tStart\tEnd\n")
        preChr = ""
        tStart = 1
        tEnd = 1
        for line in nf:
            if line.startswith('#'):
                continue
            line = line.strip()
            arr = line.split("\t")
            if arr[5] != "0":
                #tchr = preChr.split(sep)[-1]
                chrArr = preChr.split(sep)
                tchr = chrArr[-1]
                cl.write(tchr + "\t" + tStart + "\t" + tEnd + "\n")
                refAsm = chrArr[0] + sep + chrArr[1]
                return refAsm
                #break
            else:
                if preChr == "":
                    tStart = arr[2]
                    tEnd = arr[3]
                else:
                    if arr[1] != preChr:
                        tchr = preChr.split(sep)[-1]
                        cl.write(tchr + "\t" + tStart + "\t" + tEnd + "\n")
                        tStart = arr[2]
                        tEnd = arr[3]
                    else:
                        tEnd = arr[3]
            preChr = arr[1]

def simpPath(nodeFile,gafList,minMQ,asmLFile,pathDir):
    with open(nodeFile) as nf:
        nodeSize = {}
        for line in nf:
            if line.startswith('#'):
                continue
            line = line.strip()
            arr = line.split("\t")
            nodeSize[arr[0]] = int(arr[4])
    
    #asmSet = set()
    covInfo = {}
    pattern = re.compile("[><]")
    pat = re.compile("\s+")
    pat2 = re.compile("^\s*$")
    nAsm = 0
    with open(gafList) as gl,open(asmLFile,'w') as af:
        for line in gl:
            if line.startswith('#'):
                continue
            if pat2.match(line):
                continue
            line = line.strip()
            arr = pat.split(line)
            covInfo[arr[0]] = {}
            af.write(arr[0] + "\n")
            pathFile = os.path.join(pathDir,str(nAsm) + ".path")
            pNameFile = os.path.join(pathDir,str(nAsm) + ".name")
            nAsm += 1
            with open(arr[1]) as gf,open(pathFile,'w') as pf,open(pNameFile,'w') as nf:
                pathDict = {}
                for mapinfo in gf:
                    mapinfo = mapinfo.strip()
                    mapArr = mapinfo.split("\t")
                    if int(mapArr[11]) > minMQ:
                        mapStr = mapArr[5].replace("s","")
                        #pf.write(mapArr[0]+"\t"+mapStr+"\n")
                        if mapArr[0] not in pathDict:
                            pathDict[mapArr[0]] = {}
                        pathDict[mapArr[0]][int(mapArr[2])] = mapStr + "\t" + mapArr[2] + "\t" + mapArr[7] + "\t" + mapArr[18].split(":")[2]
                        
                for k in pathDict:
                    sortPath = sorted(pathDict[k].items(),key=lambda x : x[0])
                    for tpath in sortPath:
                        pf.write(k+"\t"+tpath[1]+"\n")
                        nf.write(k+"\n")
            
def nodeCov(nodeFile,gafList,minMQ,covFile,asmLFile,pathDir):
    with open(nodeFile) as nf:
        nodeSize = {}
        for line in nf:
            if line.startswith('#'):
                continue
            line = line.strip()
            arr = line.split("\t")
            nodeSize[arr[0]] = int(arr[4])
    
    #asmSet = set()
    covInfo = {}
    pattern = re.compile("[><]")
    pat = re.compile("\s+")
    pat2 = re.compile("^\s*$")
    nAsm = 0
    with open(gafList) as gl,open(asmLFile,'w') as af:
        for line in gl:
            if line.startswith('#'):
                continue
            if pat2.match(line):
                continue
            line = line.strip()
            arr = pat.split(line)
            covInfo[arr[0]] = {}
            af.write(arr[0] + "\n")
            pathFile = os.path.join(pathDir,str(nAsm) + ".path")
            pNameFile = os.path.join(pathDir,str(nAsm) + ".name")
            nAsm += 1
            with open(arr[1]) as gf,open(pathFile,'w') as pf,open(pNameFile,'w') as nf:
                pathDict = {}
                for mapinfo in gf:
                    mapinfo = mapinfo.strip()
                    mapArr = mapinfo.split("\t")
                    if int(mapArr[11]) > minMQ:
                        mapStr = mapArr[5].replace("s","")
                        #pf.write(mapArr[0]+"\t"+mapStr+"\n")
                        if mapArr[0] not in pathDict:
                            pathDict[mapArr[0]] = {}
                        pathDict[mapArr[0]][int(mapArr[2])] = mapStr + "\t" + mapArr[2] + "\t" + mapArr[7] + "\t" + mapArr[18].split(":")[2]
                        
                        nodeArr = pattern.split(mapStr)
                        nodeCount = len(nodeArr)
                        if nodeCount < 3:
                            #covInfo[arr[0]][nodeArr[1]] = (int(mapArr[8]) - int(mapArr[7]) + 1) / nodeSize[nodeArr[1]]
                            covInfo[arr[0]][nodeArr[1]] = (int(mapArr[8]) - int(mapArr[7])) / nodeSize[nodeArr[1]]
                        else:
                            for i in range(1,nodeCount):
                                if i == 1:
                                    if nodeArr[1] in covInfo[arr[0]]:
                                        covInfo[arr[0]][nodeArr[1]] += (nodeSize[nodeArr[1]] - int(mapArr[7])) / nodeSize[nodeArr[1]]
                                    else:
                                        covInfo[arr[0]][nodeArr[1]] = (nodeSize[nodeArr[1]] - int(mapArr[7])) / nodeSize[nodeArr[1]]
                                elif i == nodeCount - 1:
                                    if nodeArr[i] in covInfo[arr[0]]:
                                        #covInfo[arr[0]][nodeArr[i]] += (nodeSize[nodeArr[i]] + int(mapArr[8]) + 1 - int(mapArr[6])) / nodeSize[nodeArr[i]]
                                        covInfo[arr[0]][nodeArr[i]] += (nodeSize[nodeArr[i]] + int(mapArr[8]) - int(mapArr[6])) / nodeSize[nodeArr[i]]
                                    else:
                                        #covInfo[arr[0]][nodeArr[i]] = (nodeSize[nodeArr[i]] + int(mapArr[8]) + 1 - int(mapArr[6])) / nodeSize[nodeArr[i]]
                                        covInfo[arr[0]][nodeArr[i]] = (nodeSize[nodeArr[i]] + int(mapArr[8]) - int(mapArr[6])) / nodeSize[nodeArr[i]]

                                else:
                                    if nodeArr[i] in covInfo[arr[0]]:
                                        covInfo[arr[0]][nodeArr[i]] += 1.00
                                    else:
                                        covInfo[arr[0]][nodeArr[i]] = 1.00
                for k in pathDict:
                    sortPath = sorted(pathDict[k].items(),key=lambda x : x[0])
                    for tpath in sortPath:
                        pf.write(k+"\t"+tpath[1]+"\n")
                        nf.write(k+"\n")
                    
    allNodes = nodeSize.keys()
    allAsm = covInfo.keys()
    numLimit = len(allAsm) / 2
    with open(covFile,"w") as cf:
        #strAsm = "\t".join(allAsm)
        #cf.write("#Segid\t" + strAsm + "\n")
        for tnode in allNodes:
            oneList = []
            thList = []
            vthList = []
            i = 0
            allValue = []
            for asmb in allAsm:
                if tnode in covInfo[asmb]:
                    if covInfo[asmb][tnode] > 0.99 and covInfo[asmb][tnode] < 1.01:
                        oneList.append(i)
                    elif covInfo[asmb][tnode] > 0:
                        thList.append(i)
                        vthList.append("%.2f" % covInfo[asmb][tnode])
                    allValue.append("%.2f" % covInfo[asmb][tnode])
                else:
                    allValue.append(0.00)
                i += 1
            
            strCov = ""
            num = len(oneList) + len(thList)
            if num < numLimit:
                if len(oneList) > 0:
                    strCov += ",".join([str(x) for x in oneList])
                else:
                    strCov = "*"
                if len(thList) > 0:
                    strCov += ("\t" + ",".join([str(x) for x in thList]))
                    strCov += ("\t" + ",".join([str(x) for x in vthList]))
            else:
                strCov = ",".join([str(x) for x in allValue])
            cf.write(tnode + "\t" + strCov + "\n")

def wSep(sep,sepFile): 
    with open(sepFile,'w') as sf:
        sf.write(sep + "\n")

def fromScratch(ncalCov,sep,asmList,minigraph,graphOpt,thread,minMQ,outDir):
    asmInfo = changeHeader(sep,asmList,outDir)
    rgraph(asmInfo["dotList"],asmInfo["suffix"],minigraph,graphOpt,thread,outDir)
    mapAsm(asmInfo["nameList"],asmInfo["dotList"],asmInfo["suffix"],minigraph,thread,outDir)
    
    upDir = os.path.join(outDir,"upload")
    os.mkdir(upDir)
    
    nodeFile = os.path.join(upDir,"node.info")
    edgeFile = os.path.join(upDir,"edge.info")
    chrListFile = os.path.join(upDir,"chr.list")
    comChrFile = os.path.join(upDir,"complete.chr.list")
    asmListFile = os.path.join(upDir,"asm.list")
    ndAsmLFile = os.path.join(upDir,"node.asm.list")
    covFile = os.path.join(upDir,"cover.info")
    sepFile = os.path.join(upDir,"sep.info")
    pathDir = os.path.join(upDir,"path")
    os.mkdir(pathDir)
    
    graphDir = os.path.join(outDir,"newGraph")
    gfaFile = os.path.join(graphDir,"out.gfa")
    
    mapDir = os.path.join(outDir,"mapDir")
    mapListFile = os.path.join(mapDir,"gaf.list")
    
    wSep(sep,sepFile)
    reduceGFA(sep,gfaFile,nodeFile,edgeFile,comChrFile,ndAsmLFile)
    refAsm = refChrList(sep,nodeFile,chrListFile)
    if ncalCov:
        simpPath(nodeFile,mapListFile,minMQ,asmListFile,pathDir)
    else:
        nodeCov(nodeFile,mapListFile,minMQ,covFile,asmListFile,pathDir)
    #
    refGFA = os.path.join(upDir,"input.ref.gfa")
    mvCommand = "mv {} {}".format(gfaFile,refGFA)
    os.system(mvCommand)

def fromGFA(ncalCov,sep,gfaFile,mapListFile,minMQ,outDir):
    upDir = os.path.join(outDir,"upload")
    os.mkdir(upDir)
    
    nodeFile = os.path.join(upDir,"node.info")
    edgeFile = os.path.join(upDir,"edge.info")
    chrListFile = os.path.join(upDir,"chr.list")
    comChrFile = os.path.join(upDir,"complete.chr.list")
    asmListFile = os.path.join(upDir,"asm.list")
    ndAsmLFile = os.path.join(upDir,"node.asm.list")
    covFile = os.path.join(upDir,"cover.info")
    sepFile = os.path.join(upDir,"sep.info")
    pathDir = os.path.join(upDir,"path")
    os.mkdir(pathDir)
    
    wSep(sep,sepFile)
    reduceGFA(sep,gfaFile,nodeFile,edgeFile,comChrFile,ndAsmLFile)
    refAsm = refChrList(sep,nodeFile,chrListFile)
    if mapListFile is not None:
        if ncalCov:
            simpPath(nodeFile,mapListFile,minMQ,asmListFile,pathDir)
        else:
            nodeCov(nodeFile,mapListFile,minMQ,covFile,asmListFile,pathDir)
    #
    refGFA = os.path.join(upDir,"input.ref.gfa")
    if gfaFile.endswith(".gz"):
        cpCommand = "zcat {} > {}".format(gfaFile,refGFA)
    else:
        cpCommand = "cp {} {}".format(gfaFile,refGFA)
    os.system(cpCommand)
    
        
def indexGraph(outDir,nthread):
    upDir = os.path.join(outDir,"upload")
    mp = minipg.GraphRange(upDir,0)
    rangeSize = paras.range
    storeDep = paras.xDep
    ex = 1000000
    spChrFile = "00000000"
    mp.edgeWrite(spChrFile,rangeSize,ex,0,nthread,storeDep)    
 
def miniMain():
    sep = paras.sep
    gfaFile = paras.rGFA
    mapListFile = paras.gafList
    minMQ = paras.minMQ
    outDir = paras.outDir
    asmList = paras.asmList
    minigraph = paras.minigraph
    graphOpt = paras.graphOpt
    thread = str(paras.thread)
    nthread = paras.thread
    ncalCov = paras.ncov
    
    checkOut(outDir)
    if gfaFile is not None:
        fromGFA(ncalCov,sep,gfaFile,mapListFile,minMQ,outDir)
    else:
        if asmList is not None and minigraph is not None:
            fromScratch(ncalCov,sep,asmList,minigraph,graphOpt,thread,minMQ,outDir)
        else:
            print("Error: lack of parameters!")
            exit(1)
    #
    if paras.index:
        indexGraph(outDir,nthread)
    
if __name__ == '__main__':
    miniMain()














import argparse
import re
import os
import gzip

parser = argparse.ArgumentParser()

parser.add_argument('--sep',help="delimiter between sample and haplotype names, by default: '#'",default="#")

parser.add_argument('--rGFA',help="Input rGFA file",default=None)
parser.add_argument('--gafList',help="Input GAF file list",default=None)
parser.add_argument('--minMQ',help="Minimal mapping quality that are allowed when process GAF files, by default: 0",type=int,default=0)
parser.add_argument('--outDir',help="Output directory",default=None,required=True)

parser.add_argument('--minigraph',help="The absolute file path of minigraph executable file, this is used when rGFA file does not exist",default=None)
parser.add_argument('--assList',help="Input assembly file list. When rGFA file does not exist this option can be taken to build the reference pangenome graph",default=None)
parser.add_argument('--thread',help="Threads of minigraph running, by default: 10",type=int,default=10)
parser.add_argument('--graphOpt',help="minigraph options (except '-cxggs -t') to generate graph, e.g.: '-q 5 -l 100000 -L 50'",default=None)

paras = parser.parse_args()


def checkOut(outDir):
    if os.path.exists(outDir):
        print("Error: output directory has already existed! {}".format(outDir))
        exit(1)
    try:
        os.mkdir(outDir)
    except:
        print("Error: output directory create failed! {}".format(outDir))
                
def changeHeader(sep,assList,outDir):  
    newDir = os.path.join(outDir,"newAssDir")
    os.mkdir(newDir)
    
    pat = re.compile("\s+")
    pat2 = re.compile("^\s*$")
    nameList = []
    dotList = []
    suffix = []
    with open(assList) as al:
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
    newAssDir = os.path.join(outDir,"newAssDir")
    newAssList = [ os.path.join(newAssDir, x + y) for x,y in zip(dotList,suffix)]
    allFa = " ".join(newAssList)
    
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
        print("Error: Graph generation failed. Please check the options!")
    
def mapAss(nameList,dotList,suffix,minigraph,thread,outDir):
    newAssDir = os.path.join(outDir,"newAssDir")
    newAssList = [ os.path.join(newAssDir, x + y) for x,y in zip(dotList,suffix)]
    
    mapDir = os.path.join(outDir,"mapDir")
    os.mkdir(mapDir)
    mapListFile = os.path.join(mapDir,"gaf.list")
    
    graphDir = os.path.join(outDir,"newGraph")
    gfaFile = os.path.join(graphDir,"out.gfa")
    
    with open(mapListFile,'w') as fh:
        #for faFile,name in zip(newAssList,dotList):
        endPoint = len(dotList)
        for i in range(endPoint):
            outGAF = os.path.join(mapDir,dotList[i] + ".gaf")
            command = minigraph + " -cxasm -t " + thread + " --vc " + gfaFile + " " + newAssList[i] + " -o " + outGAF
            os.system(command)
            
            fh.write(nameList[i] + "\t" + outGAF + "\n")

def reduceGFA(gfaFile,nodeFile,edgeFile):
    flag = False
    if gfaFile.endswith(".gz"):
        fh = gzip.open(gfaFile,"rb")
        flag = True
    else:
        fh = open(gfaFile)
    with open(nodeFile,'w') as nfh, open(edgeFile,'w') as efh:
        nfh.write("#Segment\tChr\tStart\tEnd\tLen\tRefOrNot\n")
        efh.write("#Source\tTarget\tOrigin1\tOrigin2\n")
        for line in fh:
            if flag:
                line = bytes.decode(line)
            line = line.strip()
            arr = line.split("\t")
            if arr[0] == 'S':
                segment = arr[1]
                slen = arr[3].split(":")[2]
                schr = arr[4].split(":")[2]
                start = arr[5].split(":")[2]
                off_start = str(int(start) + 1)
                end = str(int(start) + int(slen))
                refOrNot = arr[6].split(":")[2]
                nfh.write(segment + "\t" + schr + "\t" + off_start + "\t" + end +"\t" + slen + "\t" + refOrNot +"\n")
            elif arr[0] == 'L':
                source = arr[1]
                origin1 = arr[2]
                target = arr[3]
                origin2 = arr[4]
                efh.write(source + "\t" + target + "\t" + origin1 + "\t" + origin2 + "\n")
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
                tchr = preChr.split(sep)[-1]
                cl.write(tchr + "\t" + tStart + "\t" + tEnd + "\n")
                break
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
    
def nodeCov(nodeFile,gafList,minMQ,covFile,assLFile,pathFile):
    with open(nodeFile) as nf:
        nodeSize = {}
        for line in nf:
            if line.startswith('#'):
                continue
            line = line.strip()
            arr = line.split("\t")
            nodeSize[arr[0]] = int(arr[4])
    
    covInfo = {}
    pattern = re.compile("[><]")
    pat = re.compile("\s+")
    pat2 = re.compile("^\s*$")
    with open(gafList) as gl,open(assLFile,'w') as al,open(pathFile,'w') as pf:
        for line in gl:
            if line.startswith('#'):
                continue
            if pat2.match(line):
                continue
            line = line.strip()
            arr = pat.split(line)
            covInfo[arr[0]] = {}
            al.write(arr[0]+"\n")
            
            with open(arr[1]) as gf:
                for mapinfo in gf:
                    mapinfo = mapinfo.strip()
                    mapArr = mapinfo.split("\t")
                    if int(mapArr[11]) > minMQ:
                        pf.write(mapArr[0]+"\t"+mapArr[5]+"\n")
                        
                        nodeArr = pattern.split(mapArr[5])
                        nodeCount = len(nodeArr)
                        if nodeCount < 3:
                            covInfo[arr[0]][nodeArr[1]] = (int(mapArr[8]) - int(mapArr[7]) + 1) / nodeSize[nodeArr[1]]
                        else:
                            for i in range(1,nodeCount):
                                if i == 1:
                                    if nodeArr[1] in covInfo[arr[0]]:
                                        covInfo[arr[0]][nodeArr[1]] += (nodeSize[nodeArr[1]] - int(mapArr[7])) / nodeSize[nodeArr[1]]
                                    else:
                                        covInfo[arr[0]][nodeArr[1]] = (nodeSize[nodeArr[1]] - int(mapArr[7])) / nodeSize[nodeArr[1]]
                                elif i == nodeCount - 1:
                                    if nodeArr[i] in covInfo[arr[0]]:
                                        covInfo[arr[0]][nodeArr[i]] += (nodeSize[nodeArr[i]] + int(mapArr[8]) + 1 - int(mapArr[6])) / nodeSize[nodeArr[i]]
                                    else:
                                        covInfo[arr[0]][nodeArr[i]] = (nodeSize[nodeArr[i]] + int(mapArr[8]) + 1 - int(mapArr[6])) / nodeSize[nodeArr[i]]

                                else:
                                    if nodeArr[i] in covInfo[arr[0]]:
                                        covInfo[arr[0]][nodeArr[i]] += 1.00
                                    else:
                                        covInfo[arr[0]][nodeArr[i]] = 1.00
    
    allNodes = nodeSize.keys()
    allAss = covInfo.keys()
    with open(covFile,"w") as cf:
        strAss = "\t".join(allAss)
        cf.write("#Segid\t" + strAss + "\n")
        for tnode in allNodes:
            covList = []
            for ass in allAss:
                if tnode in covInfo[ass]:
                    covList.append("%.2f" % covInfo[ass][tnode])
                else:
                    covList.append("0.00")
            strCov = "\t".join(covList)
            cf.write(tnode + "\t" + strCov + "\n")    

def wSep(sep,sepFile): 
    with open(sepFile,'w') as sf:
        sf.write(sep + "\n")

def fromScratch(sep,assList,minigraph,graphOpt,thread,minMQ,outDir):
    assInfo = changeHeader(sep,assList,outDir)
    rgraph(assInfo["dotList"],assInfo["suffix"],minigraph,graphOpt,thread,outDir)
    mapAss(assInfo["nameList"],assInfo["dotList"],assInfo["suffix"],minigraph,thread,outDir)
    
    upDir = os.path.join(outDir,"upload")
    os.mkdir(upDir)
    
    nodeFile = os.path.join(upDir,"node.info")
    edgeFile = os.path.join(upDir,"edge.info")
    chrListFile = os.path.join(upDir,"chr.list")
    assListFile = os.path.join(upDir,"ass.list")
    covFile = os.path.join(upDir,"cover.info")
    pathFile = os.path.join(upDir,"path.info")
    sepFile = os.path.join(upDir,"sep.info")
    
    graphDir = os.path.join(outDir,"newGraph")
    gfaFile = os.path.join(graphDir,"out.gfa")
    
    mapDir = os.path.join(outDir,"mapDir")
    mapListFile = os.path.join(mapDir,"gaf.list")
    
    wSep(sep,sepFile)
    reduceGFA(gfaFile,nodeFile,edgeFile)
    refChrList(sep,nodeFile,chrListFile)
    nodeCov(nodeFile,mapListFile,minMQ,covFile,assListFile,pathFile)

def fromGFA(sep,gfaFile,mapListFile,minMQ,outDir):
    upDir = os.path.join(outDir,"upload")
    os.mkdir(upDir)
    
    nodeFile = os.path.join(upDir,"node.info")
    edgeFile = os.path.join(upDir,"edge.info")
    chrListFile = os.path.join(upDir,"chr.list")
    assListFile = os.path.join(upDir,"ass.list")
    covFile = os.path.join(upDir,"cover.info")
    pathFile = os.path.join(upDir,"path.info")
    sepFile = os.path.join(upDir,"sep.info")
    
    wSep(sep,sepFile)
    reduceGFA(gfaFile,nodeFile,edgeFile)
    refChrList(sep,nodeFile,chrListFile)
    if mapListFile is not None:
        nodeCov(nodeFile,mapListFile,minMQ,covFile,assListFile,pathFile)
 
def miniMain():
    sep = paras.sep
    gfaFile = paras.rGFA
    mapListFile = paras.gafList
    minMQ = paras.minMQ
    outDir = paras.outDir
    assList = paras.assList
    minigraph = paras.minigraph
    graphOpt = paras.graphOpt
    thread = str(paras.thread)
    checkOut(outDir)
    if gfaFile is not None:
        fromGFA(sep,gfaFile,mapListFile,minMQ,outDir)
    else:
        if assList is not None and minigraph is not None:
            fromScratch(sep,assList,minigraph,graphOpt,thread,minMQ,outDir)
        else:
            if assList is not None and os.access(minigraph,os.X_OK):
                fromScratch(sep,assList,minigraph,graphOpt,thread,minMQ,outDir)
            else:
                print("Error: lack of parameters!")
                exit(1)
    
if __name__ == '__main__':
    miniMain()












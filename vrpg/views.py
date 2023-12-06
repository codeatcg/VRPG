from django.shortcuts import render
from django.http import JsonResponse
from django.views.decorators.csrf import csrf_exempt


import os
import re
#
from module import minipg

BinDir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
def index(request):
    allFiles = os.listdir(os.path.join(BinDir,"upload"))
    allDir = [i  for i in allFiles if os.path.isdir(os.path.join(BinDir,"upload",i)) and i != "path" and i != "anno"]
    return render(request,"vrpg/index.html",{"folder":allDir})

def showInfo(request):
    return render(request,"vrpg/info.html",{})
    
def showManual(request):
    return render(request,"vrpg/manual.html",{})       

def getSep(sepFile):
    with open(sepFile) as sf:
        line = sf.read().strip()
    return line

def readChr(chrFile):
    nameList = []
    lenList = []
    with open(chrFile) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            line = line.strip()
            arr = line.split("\t")
            nameList.append(arr[0])
            lenList.append(arr[2])
    return {"nameList":nameList,"lenList":lenList}        

def readAss(assFile):
    assList = []
    with open(assFile) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            line = line.strip()
            assList.append(line)
    return assList

def qNodesCov(covFile,nodeVec,covNameFile,ass):
    assPos = 0
    pos = 0
    with open(covNameFile) as af:
        for line in af:
            line = line.strip()
            if line == ass:
                assPos = pos
            pos += 1
    
    nodeDict = {}
    nodeSort = sorted(set(nodeVec))
    nodePos = 0
    num = len(nodeSort)
    cov = []
    flag = False
    with open(covFile) as cf:
        for line in cf:
            line = line.strip()
            arr = line.split("\t")
            if nodeSort[nodePos] == int(arr[0]):
                flag = False
                if arr[1] != "*":
                    oneArr = arr[1].split(",")
                    if len(oneArr) < pos:
                        for i in oneArr:
                            if int(i) == assPos:
                                nodeDict[nodeSort[nodePos]] = 1.00
                                flag = True
                                break
                    else:
                        nodeDict[nodeSort[nodePos]] = oneArr[assPos]
                        flag = True
                
                if not flag:
                    if len(arr) > 2:
                        mDx = arr[2].split(",")
                        mArr = arr[3].split(",")
                        for k,x in zip(mDx,mArr):
                            if int(k) == assPos:
                                nodeDict[nodeSort[nodePos]] = x
                                flag = True
                                break
                if not flag:
                    print("Error: node {} may be not in the assembly!".format(nodeSort[nodePos]))
                    nodeDict[nodeSort[nodePos]] = 0.00
                nodePos += 1
            if nodePos == num:
                break
                
    for nd in nodeVec:
        cov.append(nodeDict[nd])
    return cov

@csrf_exempt    
def showGraph(request):
    
    para = request.POST
    species = para.get("species")
    upDir = ""
    if species == '0':
        upDir = os.path.join(BinDir,"upload")
    else:
        upDir = os.path.join(BinDir,"upload",species)

    bEdgeFile = os.path.join(upDir,"edge.bw")
    eIndexFile = os.path.join(upDir,"edge.bdx")
    rNdDxFile = os.path.join(upDir,"node.ref.bdx")
    rNdFile = os.path.join(upDir,"node.ref.bw")
    nrNdFile = os.path.join(upDir,"node.nonref.bw")
    mNdDxFile = os.path.join(upDir,"node.merge.bdx")

    formFile = os.path.join(upDir,"form.info")
    chrListFile = os.path.join(upDir,"load.chr.list")
    assListFile = os.path.join(upDir,"ass.list")
    
    
    sChr = para.get("tchr")
    sStart = int(para.get("start"))
    sEnd = int(para.get("end"))
    ex = 1000000
    wStart = 50
    wWidth = 800
    wCut = 2000
    y = 300
    
    ass = para.get("ass")
    chrList = {}
    if os.path.exists(chrListFile):
        chrList = readChr(chrListFile)
    else:
        chrListFile = os.path.join(upDir,"chr.list")
        chrList = readChr(chrListFile)
        
    assList = readAss(assListFile)
    indexFlag = 0
    if os.path.exists(bEdgeFile) and os.path.exists(eIndexFile) and os.path.exists(rNdDxFile) and os.path.exists(rNdFile) and os.path.exists(nrNdFile) and os.path.exists(mNdDxFile):
        indexFlag = 1
    
    wsim = para.get("sim")
    sim = False
    refSim = False
    if wsim == "mnr":
        sim = True
        refSim = True
    elif wsim == "mr":
        sim = True
    
    queryDep = int(para.get("shdep"))
    
    mp = minipg.GraphRange(upDir,indexFlag)
    mp.formatGraph(ass,sChr,sStart,sEnd,ex,wStart,wWidth,wCut,y,queryDep,sim,refSim)
    
    draw_node = mp.draw_node
    draw_pos = mp.draw_pos
    
    ndGenePos = mp.ndGenePos
    geneVec = mp.geneVec
    layerVec = mp.layerVec
    strandVec = mp.strandVec
    
    rNodeNum = len(draw_pos)
    
    layout = para.get("lay")
    if len(draw_node) > 1000:
        layout = "ex"
    
    if layout == "co":
        for i in range(rNodeNum ):
            draw_node[i]["x"] = draw_pos[i];
            draw_node[i]["y"] = y;
            draw_node[i]["fixed"] = '1';
    else:
        for i in range(rNodeNum ):
            draw_node[i]["fx"] = draw_pos[i];
            draw_node[i]["fy"] = y;
    
    draw_edge = mp.draw_edge
    neStart = rNodeNum - 1
    dnode_len = mp.dnode_len
    for k in dnode_len:
        draw_edge[neStart]["dis"] = k
        neStart += 1
    #
    hnCov = []
    if ass != "0":
        if len(mp.hnGroup) > 0:
            nodeVec = [mp.nnames[x] for x in mp.hnGroup]
            covFile = os.path.join(upDir,"cover.info")
            covNameFile = os.path.join(upDir,"ass.list")
            if os.path.exists(covFile):
                hnCov = qNodesCov(covFile,nodeVec,covNameFile,ass)
            else:
                mm = minipg.QueryNode(upDir)
                mm.queryAssCov(nodeVec,ass)
                hnCov = mm.ndCov
            
    graphInfo = {'nodes':draw_node,'links':draw_edge, 'genome':mp.genome,'nnames':mp.nnames,'hnGroup':mp.hnGroup,'hLinks':mp.hLinks,'hDir':mp.hDir,'hnCov':hnCov,'nameList':chrList['nameList'],'lenList':chrList['lenList'],'ass':assList,
        'genePos':mp.ndGenePos,'geneVec':mp.geneVec,'layerVec':mp.layerVec,'strand':mp.strandVec
    }
    
    return JsonResponse(graphInfo)

@csrf_exempt
def initGraph(request):
    para = request.POST
    species = para.get("species")
    upDir = ""
    if species == '0':
        upDir = os.path.join(BinDir,"upload")
    else:
        upDir = os.path.join(BinDir,"upload",species)
    
    bEdgeFile = os.path.join(upDir,"edge.bw")
    eIndexFile = os.path.join(upDir,"edge.bdx")
    rNdDxFile = os.path.join(upDir,"node.ref.bdx")
    rNdFile = os.path.join(upDir,"node.ref.bw")
    nrNdFile = os.path.join(upDir,"node.nonref.bw")
    mNdDxFile = os.path.join(upDir,"node.merge.bdx")

    formFile = os.path.join(upDir,"form.info")
    chrListFile = os.path.join(upDir,"load.chr.list")
    assListFile = os.path.join(upDir,"ass.list")
    
    chrList = {}
    if os.path.exists(chrListFile):
        chrList = readChr(chrListFile)
    else:
        chrListFile = os.path.join(upDir,"chr.list")
        chrList = readChr(chrListFile)
    
    assList = readAss(assListFile)
    
    sChr = chrList["nameList"][0]
    sStart = 1
    sEnd = 10000
    if os.path.exists(formFile):
        sEnd = 1000
    
    ex = 1000000
    wStart = 50
    wWidth = 800
    wCut = 2000
    y = 300
    
    ass = "0"
    indexFlag = 0
    if os.path.exists(bEdgeFile) and os.path.exists(eIndexFile) and os.path.exists(rNdDxFile) and os.path.exists(rNdFile) and os.path.exists(nrNdFile) and os.path.exists(mNdDxFile):
        indexFlag = 1
    
    sim = True
    refSim = False
    
    queryDep = 10
    
    mp = minipg.GraphRange(upDir,indexFlag)
    mp.formatGraph(ass,sChr,sStart,sEnd,ex,wStart,wWidth,wCut,y,queryDep,sim,refSim)
    draw_node = mp.draw_node
    draw_pos = mp.draw_pos
    
    ndGenePos = mp.ndGenePos
    geneVec = mp.geneVec
    layerVec = mp.layerVec
    strandVec = mp.strandVec
    
    rNodeNum = len(draw_pos)
    
    for i in range(rNodeNum ):
        draw_node[i]["fx"] = draw_pos[i];
        draw_node[i]["fy"] = y;
            
    draw_edge = mp.draw_edge
    neStart = rNodeNum - 1
    dnode_len = mp.dnode_len
    for k in dnode_len:
        draw_edge[neStart]["dis"] = k
        neStart += 1
    ######################
    
    hnCov = []
    
    graphInfo = {'nodes':draw_node,'links':draw_edge,'genome':mp.genome,'nnames':mp.nnames,'hnGroup':mp.hnGroup,'hLinks':mp.hLinks,'hDir':mp.hDir,'hnCov':hnCov,'nameList':chrList['nameList'],'lenList':chrList['lenList'],'ass':assList,'iniEnd':sEnd,
        'genePos':mp.ndGenePos,'geneVec':mp.geneVec,'layerVec':mp.layerVec,'strand':mp.strandVec
    }
    
    return JsonResponse(graphInfo)    
    
@csrf_exempt            
def searchNode(request):
    para = request.POST
    species = para.get("species")
    node = para.get('seg')
    dbDir = ""
    if species == '0':
        dbDir = os.path.join(BinDir,"upload")
    else:
        dbDir = os.path.join(BinDir,"upload",species)
    
    dbNodeFile = os.path.join(dbDir,"node.sort.bw")
    dbCovFile = os.path.join(dbDir,"cover.bw")
    assListFile = os.path.join(dbDir,"ass.list")
    
    sepFile = os.path.join(dbDir,"sep.info")
    
    nodeAss = ''
    nodeChr = ''
    nodeStart = ''
    nodeEnd = ''
    header = []
    sep = getSep(sepFile)
    
    with open(assListFile) as af:
        for line in af:
            header.append(line.strip())
    
    geneList = [[]]
    if os.path.exists(dbNodeFile):
        mp = minipg.QueryNode(dbDir)
        mp.queryDbNode(int(node))
        
        nodeAss = mp.nodeAss
        nodeChr = mp.nodeChr
        nodeStart = mp.nodeStart
        nodeEnd = mp.nodeEnd
        
        mp.queryGene(int(node),nodeAss)
        geneList = mp.geneList
    else:
        nodeFile = os.path.join(dbDir,"node.info")
        with open(nodeFile) as nf:        
            for line in nf:
                if line.startswith('#'):
                    continue
                line = line.strip()
                arr = line.split('\t')
                if node == arr[0]:
                    nodeAssArr = arr[1].split(sep)
                    nodeAss = nodeAssArr[0] + sep + nodeAssArr[1]
                    nodeChr = nodeAssArr[-1]
                    nodeStart = arr[2]
                    nodeEnd = arr[3]
                    break
        mp = minipg.QueryNode(dbDir)
        mp.queryGene(int(node),nodeAss)
        geneList = mp.geneList
    
    
    cov = []    
    if os.path.exists(dbCovFile):
        mp = minipg.QueryNode(dbDir)
        mp.queryDbCov(int(node))
        cov = mp.ndCov
    else:
        covFile = os.path.join(dbDir,"cover.info")
        cov = [0 for i in header]
        with open(covFile) as cf:
            for line in cf:
                line = line.strip()
                arr = line.split("\t")
                if node == arr[0]:
                    if arr[1] != "*":
                        oneArr = arr[1].split(",")
                        if len(oneArr) < len(header):
                            for i in oneArr:
                                cov[int(i)] = 1.00;
                        else:
                            cov = oneArr
                            break
                            
                    if len(arr) > 2:
                        mDx = arr[2].split(",")
                        mArr = arr[3].split(",")
                        for k,x in zip(mDx,mArr):
                            cov[int(k)] = x
                    break
        
        
    return JsonResponse({'ass':header,'cov':cov,'nodeAss':nodeAss,'nodeChr':nodeChr,'nodeStart':nodeStart,'nodeEnd':nodeEnd,'geneList':geneList})   



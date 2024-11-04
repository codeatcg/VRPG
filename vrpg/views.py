from django.shortcuts import render
from django.http import JsonResponse
from django.views.decorators.csrf import csrf_exempt

import os
import re
import time
#
from module import minipg

BinDir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
#
BinDir = BinDir
def index(request):
    allFiles = os.listdir(os.path.join(BinDir,"upload"))
    dDir = {"path","anno","mapping"}
    allDir = [i  for i in allFiles if os.path.isdir(os.path.join(BinDir,"upload",i)) and i not in dDir]
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

def readAsm(asmFile):
    asmList = []
    with open(asmFile) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            line = line.strip()
            asmList.append(line)
    return asmList

def qNodesCov(covFile,nodeVec,covNameFile,asm):
    asmPos = 0
    pos = 0
    with open(covNameFile) as af:
        for line in af:
            line = line.strip()
            if line == asm:
                asmPos = pos
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
                            if int(i) == asmPos:
                                nodeDict[nodeSort[nodePos]] = 1.00
                                flag = True
                                break
                    else:
                        nodeDict[nodeSort[nodePos]] = oneArr[asmPos]
                        flag = True
                
                if not flag:
                    if len(arr) > 2:
                        mDx = arr[2].split(",")
                        mArr = arr[3].split(",")
                        for k,x in zip(mDx,mArr):
                            if int(k) == asmPos:
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
    #print(para)
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
    asmListFile = os.path.join(upDir,"asm.list")
    
    
    sChr = para.get("tchr")
    sStart = int(para.get("start"))
    sEnd = int(para.get("end"))
    buFilt = int(para.get("buFilt"))
    ex = 1000000
    wStart = 50
    wWidth = 800
    wCut = 2000
    y = 250
    
    asm = para.get("asm")
    vseq = para.get("vseq")
    if vseq == "1":
        taskID = para.get("taskID")
        asm = "!" + taskID
    chrList = {}
    if os.path.exists(chrListFile):
        chrList = readChr(chrListFile)
    else:
        chrListFile = os.path.join(upDir,"chr.list")
        chrList = readChr(chrListFile)
        
    asmList = {}
    if os.path.exists(asmListFile):
        asmList = readAsm(asmListFile)
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
        #refSim = False
        sim = True
    queryDep = int(para.get("shdep"))
    
    '''
    depFile = os.path.join(upDir,"index.dep")
    if os.path.exists(depFile):
        with open(depFile) as dh:
            queryDep = int(dh.read().strip())
    '''
    
    mp = minipg.GraphRange(upDir,indexFlag)
    mp.formatGraph(asm,sChr,sStart,sEnd,ex,wStart,wWidth,wCut,y,queryDep,buFilt,sim,refSim)
    
    draw_node = mp.draw_node
    draw_pos = mp.draw_pos
    #
    rNodeNum = len(draw_pos)
    
    layout = para.get("lay")
        
    if layout == "cosq" or layout == "coex":
        for i in range(rNodeNum ):
            draw_node[i]["fx"] = draw_pos[i];
            draw_node[i]["fy"] = 20;
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
    if vseq != '1' and asm != "" and ',' not in asm:
        if len(mp.hnGroup) > 0:
            nodeVec = [mp.nnames[x] for x in mp.hnGroup]
            covFile = os.path.join(upDir,"cover.info")
            covNameFile = os.path.join(upDir,"asm.list")
            if os.path.exists(covFile):
                hnCov = qNodesCov(covFile,nodeVec,covNameFile,asm)
            else:
                mm = minipg.QueryNode(upDir)
                mm.queryAsmCov(nodeVec,asm)
                hnCov = mm.ndCov
       
    graphInfo = {'nodes':draw_node,'links':draw_edge, 'genome':mp.genome,'nnames':mp.nnames,'hnGroup':mp.hnGroup,'hLinks':mp.hLinks,'hDir':mp.hDir,'hnCov':hnCov,'nameList':chrList['nameList'],'lenList':chrList['lenList'],'asm':asmList,
        'genePos':mp.ndGenePos,'geneVec':mp.geneVec,'layerVec':mp.layerVec,'strand':mp.strandVec,'mgFlagVec':mp.mgFlagVec,'figScale':mp.figScale,'tickValue':mp.tickValue,'tickPos':mp.tickPos,'ndExonPos':mp.ndExonPos,'rnaVec':mp.rnaVec,
        'eLayerVec':mp.eLayerVec,'eStrandVec':mp.eStrandVec,'eNumVec':mp.eNumVec,'eFlagVec':mp.eFlagVec,'ndCDSPos':mp.ndCDSPos,'cdsVec':mp.cdsVec,'cLayerVec':mp.cLayerVec,'cNumVec':mp.cNumVec,'hEdgeAsm':mp.hEdgeAsm,'qChr':mp.qChr,'qStart':mp.qStart,'qEnd':mp.qEnd,
        'qPath':mp.qPath,'qCigar':mp.qCigar,'tkNameVec':mp.tkNameVec,'tkDesVec':mp.tkDesVec,'tkColVec':mp.tkColVec,'tkCumVec':mp.tkCumVec,'tkItem':mp.tkItem,'rBedPos':mp.rBedPos,'rBedName':mp.rBedName,'rBedLayer':mp.rBedLayer,
        'rBedScore':mp.rBedScore,'rBedStrand':mp.rBedStrand
    }
    #print(mp.rBedScore)
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
    asmListFile = os.path.join(upDir,"asm.list")
    
    chrList = {}
    if os.path.exists(chrListFile):
        chrList = readChr(chrListFile)
    else:
        chrListFile = os.path.join(upDir,"chr.list")
        chrList = readChr(chrListFile)
    
    asmList = {}
    if os.path.exists(asmListFile):
        asmList = readAsm(asmListFile)
            
    sChr = chrList["nameList"][0]
    sStart = 1
    sEnd = 10000
    if os.path.exists(formFile):
        sEnd = 1000
    
    buFilt = 50
    ex = 1000000
    wStart = 50
    wWidth = 800
    wCut = 2000
    y = 250
    
    asm = ""
    indexFlag = 0
    if os.path.exists(bEdgeFile) and os.path.exists(eIndexFile) and os.path.exists(rNdDxFile) and os.path.exists(rNdFile) and os.path.exists(nrNdFile) and os.path.exists(mNdDxFile):
        indexFlag = 1

    sim = True
    refSim = False
    
    queryDep = 10
    '''
    depFile = os.path.join(upDir,"index.dep")
    if os.path.exists(depFile):
        with open(depFile) as dh:
            queryDep = int(dh.read().strip())
    '''
    mp = minipg.GraphRange(upDir,indexFlag)
    mp.formatGraph(asm,sChr,sStart,sEnd,ex,wStart,wWidth,wCut,y,queryDep,buFilt,sim,refSim)    
    draw_node = mp.draw_node
    draw_pos = mp.draw_pos
    #
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
    
    graphInfo = {'nodes':draw_node,'links':draw_edge,'genome':mp.genome,'nnames':mp.nnames,'hnGroup':mp.hnGroup,'hLinks':mp.hLinks,'hDir':mp.hDir,'hnCov':hnCov,'nameList':chrList['nameList'],'lenList':chrList['lenList'],'asm':asmList,
        'genePos':mp.ndGenePos,'geneVec':mp.geneVec,'layerVec':mp.layerVec,'strand':mp.strandVec,'mgFlagVec':mp.mgFlagVec,'figScale':mp.figScale,'tickValue':mp.tickValue,'tickPos':mp.tickPos,'ndExonPos':mp.ndExonPos,'rnaVec':mp.rnaVec,
        'eLayerVec':mp.eLayerVec,'eStrandVec':mp.eStrandVec,'eNumVec':mp.eNumVec,'eFlagVec':mp.eFlagVec,'ndCDSPos':mp.ndCDSPos,'cdsVec':mp.cdsVec,'cLayerVec':mp.cLayerVec,'cNumVec':mp.cNumVec,'tkNameVec':mp.tkNameVec,'tkDesVec':mp.tkDesVec,
        'tkColVec':mp.tkColVec,'tkCumVec':mp.tkCumVec,'tkItem':mp.tkItem,'rBedPos':mp.rBedPos,'rBedName':mp.rBedName,'rBedLayer':mp.rBedLayer,'rBedScore':mp.rBedScore,'rBedStrand':mp.rBedStrand
    }
    #
    return JsonResponse(graphInfo)    

@csrf_exempt   
def nodeGene(request):
    para = request.POST
    species = para.get("species")
    node = para.get('seg')
    dbDir = ""
    if species == '0':
        dbDir = os.path.join(BinDir,"upload")
    else:
        dbDir = os.path.join(BinDir,"upload",species)
    dbNodeFile = os.path.join(dbDir,"node.sort.bw")
    
    sepFile = os.path.join(dbDir,"sep.info")
    
    nodeAsm = ''
    nodeChr = ''
    nodeStart = ''
    nodeEnd = ''
    sep = getSep(sepFile)
    
    geneList = []
    mp = minipg.QueryNode(dbDir)
    if os.path.exists(dbNodeFile):
        mp.queryDbNode(int(node))
        nodeAsm = mp.nodeAsm
        nodeChr = mp.nodeChr
        nodeStart = mp.nodeStart
        nodeEnd = mp.nodeEnd
    else:
        nodeFile = os.path.join(dbDir,"node.info")
        with open(nodeFile) as nf:        
            for line in nf:
                if line.startswith('#'):
                    continue
                line = line.strip()
                arr = line.split('\t')
                if node == arr[0]:
                    nodeAsmArr = arr[1].split(sep)
                    nodeAsm = nodeAsmArr[0] + sep + nodeAsmArr[1]
                    nodeChr = nodeAsmArr[-1]
                    nodeStart = arr[2]
                    nodeEnd = arr[3]
                    #print(arr)
                    break
    
    mp.queryGene(int(node),nodeAsm)
    annoList = mp.geneList
    for oneAnn in annoList:
        if oneAnn[1].endswith("gene"):
            geneList.append(oneAnn)
                
    return JsonResponse({'nodeAsm':nodeAsm,'nodeChr':nodeChr,'nodeStart':nodeStart,'nodeEnd':nodeEnd,'geneList':geneList})
    
@csrf_exempt            
def searchNode(request):
    para = request.POST
    species = para.get("species")
    node = para.get('seg')
    #
    dbDir = ""
    if species == '0':
        dbDir = os.path.join(BinDir,"upload")
    else:
        dbDir = os.path.join(BinDir,"upload",species)
    dbNodeFile = os.path.join(dbDir,"node.sort.bw")
    dbCovFile = os.path.join(dbDir,"cover.bw")
    asmListFile = os.path.join(dbDir,"asm.list")
    
    #
    sepFile = os.path.join(dbDir,"sep.info")
    
    nodeAsm = ''
    nodeChr = ''
    nodeStart = ''
    nodeEnd = ''
    nodeSeq = ''
    header = []
    sep = getSep(sepFile)
    
    #
    with open(asmListFile) as af:
        for line in af:
            header.append(line.strip())
    
    geneList = [[]]
    if os.path.exists(dbNodeFile):
        mp = minipg.QueryNode(dbDir)
        mp.queryDbNode(int(node))

        nodeAsm = mp.nodeAsm
        nodeChr = mp.nodeChr
        nodeStart = mp.nodeStart
        nodeEnd = mp.nodeEnd
        
        mp.queryGene(int(node),nodeAsm)
        geneList = mp.geneList
        #
        mp.fetchNdSeq(int(node))
        nodeSeq = mp.nodeSeq
    else:
        nodeFile = os.path.join(dbDir,"node.info")
        with open(nodeFile) as nf:        
            for line in nf:
                if line.startswith('#'):
                    continue
                line = line.strip()
                arr = line.split('\t')
                if node == arr[0]:
                    nodeAsmArr = arr[1].split(sep)
                    nodeAsm = nodeAsmArr[0] + sep + nodeAsmArr[1]
                    nodeChr = nodeAsmArr[-1]
                    nodeStart = arr[2]
                    nodeEnd = arr[3]
                    #print(arr)
                    break
        mp = minipg.QueryNode(dbDir)
        mp.queryGene(int(node),nodeAsm)
        geneList = mp.geneList
    
    
    cov = []    
    if os.path.exists(dbCovFile):
        mp = minipg.QueryNode(dbDir)
        mp.queryDbCov(int(node))
        cov = mp.ndCov
    else:
        covFile = os.path.join(dbDir,"cover.info")
        if os.path.exists(covFile):
            cov = [0 for i in header]
            with open(covFile) as cf:
                for line in cf:
                    line = line.strip()
                    arr = line.split("\t")
                    if node == arr[0]:
                        if arr[1] != "*":
                            #
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
        
        
    return JsonResponse({'asm':header,'cov':cov,'nodeAsm':nodeAsm,'nodeChr':nodeChr,'nodeStart':nodeStart,'nodeEnd':nodeEnd,'geneList':geneList,'nodeSeq':nodeSeq})   

def taskQuery(species):
    tkNum = 0
    upDir = os.path.join(BinDir,"upload")
    mDir = os.path.join(upDir,"mapping")
    if os.path.exists(mDir):
        taskList = os.listdir(mDir)
        for tdir in taskList:
            tkFile = os.path.join(mDir,tdir,"task.info")
            if not os.path.exists(tkFile):
                tkNum += 1
    
    dirList = os.listdir(upDir)
    for dir in dirList:
        if dir != "mapping":
            dirPath = os.path.join(upDir,dir)
            if os.path.isdir(dirPath):
                mDir = os.path.join(dirPath,"mapping")
                if os.path.exists(mDir):
                    taskList = os.listdir(mDir)
                    for tdir in taskList:
                        tkFile = os.path.join(mDir,tdir,"task.info")
                        if not os.path.exists(tkFile):
                            tkNum += 1
    
    return tkNum
                
def createTask(species):
    preDir = ""
    if species == "0":
        preDir = os.path.join(BinDir,"upload","mapping")
    else:
        preDir = os.path.join(BinDir,"upload",species,"mapping")
    if not os.path.exists(preDir):
        os.mkdir(preDir,0o755)
        
    taskID = "task_" + str(time.time()) + "_1"
    taskDir = os.path.join(preDir,taskID)
    num = 1
    while os.path.exists(taskDir):
        num += 1
        taskID = "task_" + str(time.time()) + "_" + str(num)
        taskDir = os.path.join(preDir,taskID)
    os.mkdir(taskDir,0o755)
    return taskID

#   
def wQuerySeq(qSeq,taskDir):
    queryFile = os.path.join(taskDir,"query.fa")
    with open(queryFile,'w') as fh:
        fh.write(qSeq)

def queryMap(minigraph,gfaFile,taskDir):
    queryFile = os.path.join(taskDir,"query.fa")
    outFile = os.path.join(taskDir,"query.gaf")
    #
    mType = " lr "
    secLine = ""
    with open(queryFile) as fh:
        firLine = fh.readline()
        if not firLine.startswith('>'):
            print("Error: query format error (not FASTA format)!")
            return
        try:
            secLine = fh.read()
        except:
            print("Error: query format error (not FASTA format)!")
            return
    
    if secLine != "":
        if len(secLine) < 200:
            
            mType = " sr "
    command = minigraph + " -cx" + mType + "--vc " + gfaFile + " " + queryFile + " -o " + outFile
    os.system(command)

def map2path(taskDir):
    gafFile = os.path.join(taskDir,"query.gaf")
    pathFile = os.path.join(taskDir,"query.path")

    with open(gafFile) as gf,open(pathFile,'w') as pf:
        pathDict = {}
        for mapinfo in gf:
            mapinfo = mapinfo.strip()
            mapArr = mapinfo.split("\t")
            #
            mapStr = mapArr[5].replace("s","")
            if mapArr[0] not in pathDict:
                pathDict[mapArr[0]] = {}
            pathDict[mapArr[0]][int(mapArr[2])] = mapStr + "\t" + mapArr[2] + "\t" + mapArr[7] + "\t" + mapArr[18].split(":")[2]
                
        for k in pathDict:
            sortPath = sorted(pathDict[k].items(),key=lambda x : x[0])
            for tpath in sortPath:
                pf.write(k+"\t"+tpath[1]+"\n")
    
    
def map2Loci(preDir,gaf2rbed,edgeFile,taskDir):
    pathFile = os.path.join(taskDir,"query.path")
    outFile = os.path.join(taskDir,"query.bed")

    rChrFile = os.path.join(preDir,"chr.list")
    rnodeFile = os.path.join(preDir,"node.ref.bw")
    dxFile = os.path.join(preDir,"node.ref.bdx")
    edgeFile = os.path.join(preDir,"edge.info")
    
    command = gaf2rbed + " --chr " + rChrFile + " --rnode " + rnodeFile + " --dxnode " + dxFile + " --edge " + edgeFile + " --path " + pathFile + " --out " + outFile
    os.system(command)
    
def readQyLoci(taskDir):
    bedFile = os.path.join(taskDir,"query.bed")
    locInfo = []
    with open(bedFile) as fh:
        for line in fh:
            line = line.strip()
            arr = line.split("\t")
            locInfo.append(arr)
    return locInfo

@csrf_exempt   
def seqQuery(request):
    para = request.POST
    querySeq = para.get("qSeq")
    species = para.get("species")
    
    upDir = os.path.join(BinDir,"upload")
    preDir = upDir
    if species != "0":
        preDir = os.path.join(upDir,species)
    
    gfaFile = ""
    edgeFile = ""
    taskDir = ""
    gfaFile = os.path.join(preDir,"input.ref.gfa")
    if not os.path.exists(gfaFile):
        return JsonResponse({"tkNum": -1})
        
    tkNum = taskQuery(species)
    if tkNum > 2:
        return JsonResponse({"tkNum": tkNum})
    taskID = createTask(species)
    
    edgeFile = os.path.join(preDir,"edge.info")
    taskDir = os.path.join(preDir,"mapping",taskID)  
    wQuerySeq(querySeq,taskDir)
    minigraph = os.path.join(BinDir,"bin","minigraph")
    
    queryMap(minigraph,gfaFile,taskDir)
    map2path(taskDir)
    gaf2rbed = os.path.join(BinDir,"module","gaf2rbed")
    #gaf2rbed = os.path.join(BinDir,"bin","gaf2rbed")
    map2Loci(preDir,gaf2rbed,edgeFile,taskDir)
    #return taskDir
    locInfo = readQyLoci(taskDir)
    #
    tkFile = os.path.join(taskDir,"task.info")
    ofh = open(tkFile,'w')
    ofh.close()
    return JsonResponse({"taskID":taskID,"locInfo":locInfo,"tkNum":1})
    







    

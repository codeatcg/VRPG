from django.shortcuts import render
from django.http import JsonResponse
from django.views.decorators.csrf import csrf_exempt


import os
import re
from module import minipg
# Create your views here.

BinDir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
def index(request):
    return render(request,"vrpg/index.html",{})


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

@csrf_exempt    
def showGraph(request):
    
    para = request.POST
    species = para.get("species")
    nodeFile = os.path.join(BinDir,"upload",species,"node.info")
    edgeFile = os.path.join(BinDir,"upload",species,"edge.info")
    chrListFile = os.path.join(BinDir,"upload",species,"chr.list")
    assListFile = os.path.join(BinDir,"upload",species,"ass.list")
    pathFile = os.path.join(BinDir,"upload",species,"path.info")
    sepFile = os.path.join(BinDir,"upload",species,"sep.info")
    
    bEdgeFile = os.path.join(BinDir,"upload",species,"edge.bw")
    eIndexFile = os.path.join(BinDir,"upload",species,"edge.dx")
    
    sChr = para.get("tchr")
    sStart = int(para.get("start"))
    sEnd = int(para.get("end"))
    ex = 1000000
    wStart = 50
    wWidth = 800
    wCut = 2000
    y = 300
    
    ass = para.get("ass")
    chrList = readChr(chrListFile)
    assList = readAss(assListFile)
    if(os.path.exists(bEdgeFile) and os.path.exists(eIndexFile)):
        mp = minipg.GraphRange(nodeFile,edgeFile,pathFile,sepFile,chrListFile,bEdgeFile,eIndexFile)
    else:
        mp = minipg.GraphRange(nodeFile,edgeFile,pathFile,sepFile)
    mp.formatGraph(ass,sChr,sStart,sEnd,ex,wStart,wWidth,wCut,y)
    
    draw_node = mp.draw_node
    draw_pos = mp.draw_pos
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
        
    graphInfo = {'nodes':draw_node,'links':draw_edge, 'genome':mp.genome,'nnames':mp.nnames,'hnGroup':mp.hnGroup,'hLinks':mp.hLinks,'nameList':chrList['nameList'],'lenList':chrList['lenList'],'ass':assList}
    return JsonResponse(graphInfo)

@csrf_exempt
def initGraph(request):
    para = request.POST
    species = para.get("species")
    nodeFile = os.path.join(BinDir,"upload",species,"node.info")
    edgeFile = os.path.join(BinDir,"upload",species,"edge.info")
    chrListFile = os.path.join(BinDir,"upload",species,"chr.list")
    assListFile = os.path.join(BinDir,"upload",species,"ass.list")
    pathFile = os.path.join(BinDir,"upload",species,"path.info")
    sepFile = os.path.join(BinDir,"upload",species,"sep.info")
    
    bEdgeFile = os.path.join(BinDir,"upload",species,"edge.bw")
    eIndexFile = os.path.join(BinDir,"upload",species,"edge.dx")
    
    chrList = readChr(chrListFile)
    assList = readAss(assListFile)
    
    sChr = chrList["nameList"][0]
    sStart = 1
    sEnd = 10000
    if species == "Homo_sapiens":
        sEnd = 50000
    ex = 1000000
    wStart = 50
    wWidth = 800
    wCut = 2000
    y = 300
    
    ass = "0"
    
    if(os.path.exists(bEdgeFile) and os.path.exists(eIndexFile)):
        mp = minipg.GraphRange(nodeFile,edgeFile,pathFile,sepFile,chrListFile,bEdgeFile,eIndexFile)
    else:
        mp = minipg.GraphRange(nodeFile,edgeFile,pathFile,sepFile)
    mp.formatGraph(ass,sChr,sStart,sEnd,ex,wStart,wWidth,wCut,y)
    
    draw_node = mp.draw_node
    draw_pos = mp.draw_pos
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
    graphInfo = {'nodes':draw_node,'links':draw_edge,'genome':mp.genome,'nnames':mp.nnames,'hnGroup':mp.hnGroup,'hLinks':mp.hLinks,'nameList':chrList['nameList'],'lenList':chrList['lenList'],'ass':assList}
    return JsonResponse(graphInfo)    
    
@csrf_exempt            
def searchNode(request):
    para = request.POST
    species = para.get("species")
    node = para.get('seg')
    covFile = os.path.join(BinDir,"upload",species,"cover.info")
    nodeFile = os.path.join(BinDir,"upload",species,"node.info")
    sepFile = os.path.join(BinDir,"upload",species,"sep.info")
    sep = getSep(sepFile)
    with open(covFile) as cf:
        header = []
        cov = []
        for line in cf:
            line = line.strip()
            arr = line.split("\t")
            if line.startswith('#'):
                header = arr[1:]
            else:
                if node == arr[0]:
                    cov = arr[1:]
                    break
    
    nodeAss = ''
    nodeChr = ''
    nodeStart = ''
    nodeEnd = ''
    with open(nodeFile) as nf:        
        for line in nf:
            if line.startswith('#'):
                continue
            line = line.strip()
            arr = line.split('\t')
            if node == arr[0]:
                nodeAssArr = arr[1].split(sep)
                nodeAss = "GRCh38.p14" if len(nodeAssArr) == 1 else nodeAssArr[0] + sep + nodeAssArr[1]
                nodeChr = nodeAssArr[-1]
                nodeStart = arr[2]
                nodeEnd = arr[3]
                break
    return JsonResponse({'ass':header,'cov':cov,'nodeAss':nodeAss,'nodeChr':nodeChr,'nodeStart':nodeStart,'nodeEnd':nodeEnd})   






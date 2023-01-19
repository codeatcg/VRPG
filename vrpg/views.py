from django.shortcuts import render
from django.http import JsonResponse
from django.views.decorators.csrf import csrf_exempt


import os
import re

# Create your views here.

BinDir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
def index(request):
    return render(request,"vrpg/index.html",{})

def getMark(ori_1,ori_2,sign='+'):
    mark = 0
    if ori_1 == sign:
        if ori_2 == sign:
            mark = 2
        else:
            mark = 3        
    else:
        if ori_2 == sign:
            mark = 4
        else:
            mark = 5
    return mark
    
def revMark(mark):
    revVal = None
    if mark == 2:
        revVal = 5
    elif mark == 5:
        revVal = 2
    else:
        revVal = mark
    return revVal
    
def parseEdge(edgeFile):
    edge_dict = {}
    with open(edgeFile) as fh:
        for line in fh:
            line = line.strip()
            arr = line.split('\t')
            mark = getMark(arr[2],arr[3])
            if arr[0] in edge_dict:
                if 'o' in edge_dict[arr[0]]:
                    edge_dict[arr[0]]['o'].append([arr[1],mark])
                else:
                    edge_dict[arr[0]]['o'] = [[arr[1],mark]]
            else:
                edge_dict[arr[0]] = {}
                edge_dict[arr[0]]['o'] = [[arr[1],mark]]
                
            if arr[1] in edge_dict:
                if 'i' in edge_dict[arr[1]]:
                    edge_dict[arr[1]]['i'].append([arr[0],mark])
                else:
                    edge_dict[arr[1]]['i'] = [[arr[0],mark]]
            else:
                edge_dict[arr[1]] = {}
                edge_dict[arr[1]]['i'] = [[arr[0],mark]]
    return edge_dict

def parseNode(nodeFile,sep,sChr,sStart,sEnd,ex):
    eStart = sStart - ex if sStart - ex > 0 else 1
    eEnd = sEnd + ex
    exNode = set()
    rangeNode = []
    info = {}
    ass = {}
    realLen = 0
    with open(nodeFile) as fh:
        fh.readline()
        for line in fh:
            line = line.strip()
            arr = line.split('\t')
            assArr = arr[1].split(sep)
            tName = "REF" if len(assArr) == 1 else assArr[0]
            if arr[5] == "0":
                tchr = assArr[-1]
                if tchr == sChr:
                    if int(arr[2]) <= eStart:
                        if int(arr[3]) >= eStart:
                            exNode.add(arr[0])
                            info[arr[0]] = int(arr[4])
                            ass[arr[0]] = tName
                            if int(arr[3]) >= sStart:
                                rangeNode.append(arr[0])
                                realLen += int(arr[4])
                    else:
                        if int(arr[2]) <= eEnd:
                            exNode.add(arr[0])
                            info[arr[0]] = int(arr[4])
                            ass[arr[0]] = tName
                        
                            if int(arr[2]) <= sStart:
                                if int(arr[3]) >= sStart:
                                    rangeNode.append(arr[0])
                                    realLen += int(arr[4])
                            else:
                                if int(arr[2]) <= sEnd:
                                    rangeNode.append(arr[0])
                                    realLen += int(arr[4])
            else:
                info[arr[0]] = int(arr[4])
                ass[arr[0]] = tName
    node_dict = {'range':rangeNode, 'ex':exNode, 'info':info, 'ass':ass, 'realLen':realLen}
    return node_dict 

def hAssNode(ass,pathFile,sep,nid_dict,r_edge_dict):
    pat1 = re.compile("[<>]")
    pat2 = re.compile("\w+")
    hnGroup = []
    hLinks = []
    #
    flag = False
    with open(pathFile) as fh:
        for line in fh:
            line = line.strip()
            arr = line.split('\t')
            tAssArr = arr[0].split(sep)
            tAss = tAssArr[0] + sep + tAssArr[1]
            if(tAss == ass):
                flag = True
                orient = pat2.split(arr[1])
                nodes = pat1.split(arr[1])
                
                nodes = nodes[1:]
                
                preOri = ""
                preIn = 0
                preNode = ""
                for k,r in zip(nodes,orient):
                    if k in nid_dict:
                        tIn = 1
                        hnGroup.append(nid_dict[k][2])
                        if preIn > 0:
                            mark = getMark(preOri,r,'>')
                            sym = preNode + "_" + k + "_" + str(mark)
                            #
                            if sym in r_edge_dict:
                                if mark == 2:
                                    if r_edge_dict[sym] == 0:
                                        mark = 1
                                if mark == 1:
                                    hLinks.append(str(nid_dict[preNode][1]) + "_" + str(nid_dict[k][0]) + "_" + str(mark))
                                elif mark == 2:
                                    hLinks.append(str(nid_dict[preNode][1]) + "_" + str(nid_dict[k][0]) + "_" + str(mark))
                                elif mark == 3:
                                    hLinks.append(str(nid_dict[preNode][1]) + "_" + str(nid_dict[k][1]) + "_" + str(mark))
                                elif mark == 4:
                                    hLinks.append(str(nid_dict[preNode][0]) + "_" + str(nid_dict[k][0]) + "_" + str(mark))
                                else:
                                    hLinks.append(str(nid_dict[preNode][0]) + "_" + str(nid_dict[k][1]) + "_" + str(mark))
                            else:
                                rmark = revMark(mark)
                                sym = k + "_" + preNode + "_" + str(rmark)
                                if sym in r_edge_dict:
                                    if rmark == 2:
                                        if r_edge_dict[sym] == 0:
                                            rmark = 1
                                    if rmark == 1:
                                        hLinks.append(str(nid_dict[k][1]) + "_" + str(nid_dict[preNode][0]) + "_" + str(rmark))
                                    elif rmark == 2:
                                        hLinks.append(str(nid_dict[k][1]) + "_" + str(nid_dict[preNode][0]) + "_" + str(rmark))
                                    elif rmark == 3:
                                        hLinks.append(str(nid_dict[k][1]) + "_" + str(nid_dict[preNode][1]) + "_" + str(rmark))
                                    elif rmark == 4:
                                        hLinks.append(str(nid_dict[k][0]) + "_" + str(nid_dict[preNode][0]) + "_" + str(rmark))
                                    else:
                                        hLinks.append(str(nid_dict[k][0]) + "_" + str(nid_dict[preNode][1]) + "_" + str(rmark))
                                
                    else:
                        tIn = 0
                    preOri = r
                    preIn = tIn
                    preNode = k
            else:
                if flag:
                    break
                
    return {"hnGroup":hnGroup,"hLinks":hLinks}                

def getSep(sepFile):
    with open(sepFile) as sf:
        line = sf.read().strip()
    return line
    
def formatGraph(nodeFile,edgeFile,pathFile,sepFile,ass,sChr,sStart,sEnd,ex,wStart,wWidth,wCut,y):
    sep = getSep(sepFile)
    edge_dict = parseEdge(edgeFile)        
    node_dict = parseNode(nodeFile,sep,sChr,sStart,sEnd,ex)
    
    realLen = node_dict['realLen']
    rcount = len(node_dict['range'])
    wPerK = wWidth / (realLen + (rcount-1) * wCut)
    x_wCut = wPerK * wCut
    
    gNum = 0
    iNum = 0
    r_edge_dict = {}
    draw_node = [{'id':0,'group':0,'fx':wStart,'fy':y}]
    draw_edge = []
    node_pre = wStart
    genome = []
    nRefNode = {}
    nid_dict = {}
    s_nid = 0
    e_nid = 0
    pre_node = ""
    range_set = set(node_dict['range'])
    nnames = []
    for node in node_dict['range']:
        node_len = node_dict['info'][node]
        if node in edge_dict:
            if 'o' in edge_dict[node]:
                tNref = []
                deep = []
                for x in edge_dict[node]['o']:
                    if x[0] in node_dict['ex']:
                        if x[0] in range_set:
                            sym = node + "_" + x[0] + "_" + str(x[1])
                            if sym not in r_edge_dict:
                                r_edge_dict[sym] = 2
                    else:
                        if x[0] not in nRefNode:
                            tNref.append(x[0])
                            deep.append(0)
                            nRefNode[x[0]] = 1
                        sym = node + "_" + x[0] + "_" + str(x[1])
                        r_edge_dict[sym] = 2
                if len(tNref) > 0:
                    i = 0
                    while i < len(tNref):
                        if deep[i] > 10:
                            break
                        if tNref[i] in edge_dict:
                            if 'o' in edge_dict[tNref[i]]:
                                for k in edge_dict[tNref[i]]['o']:
                                    if k[0] in node_dict['ex']:
                                        if k[0] in range_set:
                                            sym = tNref[i] + "_" + k[0] + "_" + str(k[1])
                                            r_edge_dict[sym] = 2
                                    else:
                                        if k[0] not in nRefNode:
                                            tNref.append(k[0])
                                            deep.append(deep[i]+1)
                                            nRefNode[k[0]] = 1
                                            
                                        sym = tNref[i] + "_" + k[0] + "_" + str(k[1])
                                        r_edge_dict[sym] = 2
                                            
                            #
                            if 'i' in edge_dict[tNref[i]]:
                                for k in edge_dict[tNref[i]]['i']:
                                    if k[0] in node_dict['ex']:
                                        if k[0] in range_set:
                                            sym = k[0] + "_" + tNref[i] + "_" + str(k[1])
                                            r_edge_dict[sym] = 2
                                    else:
                                        if k[0] not in nRefNode:
                                            tNref.append(k[0])
                                            deep.append(deep[i]+1)
                                            nRefNode[k[0]] = 1
                                            
                                        sym = k[0] + "_" + tNref[i] + "_" + str(k[1])
                                        r_edge_dict[sym] = 2
                        i += 1
            
            #
            if 'i' in edge_dict[node]:
                tNref = []
                deep = []
                for x in edge_dict[node]['i']:
                    if x[0] not in node_dict['ex']:
                        if x[0] not in nRefNode:
                            tNref.append(x[0])
                            deep.append(0)
                            nRefNode[x[0]] = 1
                        sym = x[0] + "_" + node + "_" + str(x[1])
                        r_edge_dict[sym] = 2
                if len(tNref) > 0:
                    i = 0
                    while i < len(tNref):
                        if deep[i] > 10:
                            break
                        if tNref[i] in edge_dict:
                            if 'o' in edge_dict[tNref[i]]:
                                for k in edge_dict[tNref[i]]['o']:
                                    if k[0] in node_dict['ex']:
                                        if k[0] in range_set:
                                            sym = tNref[i] + "_" + k[0] + "_" + str(k[1])
                                            r_edge_dict[sym] = 2
                                    else:
                                        if k[0] not in nRefNode:
                                            tNref.append(k[0])
                                            deep.append(deep[i]+1)
                                            nRefNode[k[0]] = 1
                                            
                                        sym = tNref[i] + "_" + k[0] + "_" + str(k[1])
                                        r_edge_dict[sym] = 2
                                            
                            #
                            if 'i' in edge_dict[tNref[i]]:
                                for k in edge_dict[tNref[i]]['i']:
                                    if k[0] in node_dict['ex']:
                                        if k[0] in range_set:
                                            sym = k[0] + "_" + tNref[i] + "_" + str(k[1])
                                            r_edge_dict[sym] = 2
                                    else:
                                        if k[0] not in nRefNode:
                                            tNref.append(k[0])
                                            deep.append(deep[i]+1)
                                            nRefNode[k[0]] = 1
                                            
                                        sym = k[0] + "_" + tNref[i] + "_" + str(k[1])
                                        r_edge_dict[sym] = 2
                        i += 1
            
        if gNum > 0:
            node_pre += x_wCut
            iNum += 1
            draw_node.append({'id':iNum,'group':gNum,'fx':node_pre,'fy':y})
            draw_edge.append({'source':(iNum-1),'target':iNum,'type':1})
            s_nid = iNum
            
            tsym = pre_node + "_" + node + "_2"
            if tsym in r_edge_dict:
                r_edge_dict[tsym] = 0
              
        if node_len <= wCut:
            trans = node_len * wPerK
            node_pre += trans
            iNum += 1
            draw_node.append({'id':iNum,'group':gNum,'fx':node_pre,'fy':y})
            draw_edge.append({'source':(iNum-1),'target':iNum,'type':0})
        else:
            preLen = wCut
            while preLen < node_len:
                node_pre += x_wCut
                iNum += 1
                draw_node.append({'id':iNum,'group':gNum,'fx':node_pre,'fy':y})
                draw_edge.append({'source':(iNum-1),'target':iNum,'type':0})
                preLen += wCut
            trans = (node_len + wCut - preLen) * wPerK
            node_pre += trans
            iNum += 1
            draw_node.append({'id':iNum,'group':gNum,'fx':node_pre,'fy':y})
            draw_edge.append({'source':(iNum-1),'target':iNum,'type':0}) 
        
        nnames.append(node)
        genome.append(node_dict['ass'][node])
        
        e_nid = iNum
        nid_dict[node] = [s_nid,e_nid,gNum]
        
        gNum += 1
        pre_node = node
    
    for nk in nRefNode:
        node_len = node_dict['info'][nk]
        if node_len <= wCut:
            iNum += 1
            s_nid = iNum
            draw_node.append({'id':iNum,'group':gNum})
            
            iNum += 1
            e_nid = iNum
            draw_node.append({'id':iNum,'group':gNum})
            draw_edge.append({'source':(iNum-1),'target':iNum,'type':0}) 
        else:
            preLen = wCut
            s_nid = iNum + 1
            while preLen < node_len:
                iNum += 1
                draw_node.append({'id':iNum,'group':gNum})
                draw_edge.append({'source':iNum,'target':(iNum+1),'type':0})
                preLen += wCut
            iNum += 1
            e_nid = iNum
            draw_node.append({'id':iNum,'group':gNum})
        
        nnames.append(nk)
        genome.append(node_dict['ass'][nk])
        
        nid_dict[nk] = [s_nid,e_nid,gNum]
        gNum += 1

    #
    for ed in r_edge_dict:
        if r_edge_dict[ed] > 1:
            t_arr = ed.split('_')
            node1 = t_arr[0]
            node2 = t_arr[1]
            if t_arr[2] == '2':

                draw_edge.append({'source':nid_dict[node1][1],'target':nid_dict[node2][0],'type':int(t_arr[2])})
            elif t_arr[2] == '3':
                draw_edge.append({'source':nid_dict[node1][1],'target':nid_dict[node2][1],'type':int(t_arr[2])})
            elif t_arr[2] == '4':
                draw_edge.append({'source':nid_dict[node1][0],'target':nid_dict[node2][0],'type':int(t_arr[2])})
            else:
                draw_edge.append({'source':nid_dict[node1][0],'target':nid_dict[node2][1],'type':int(t_arr[2])})
    
    if ass == "0":
        return {'nodes':draw_node, 'links':draw_edge, 'genome':genome,'nnames':nnames,'hnGroup':[],'hLinks':[]}
    else:
        assInfo = hAssNode(ass,pathFile,sep,nid_dict,r_edge_dict)
        #
        return {'nodes':draw_node, 'links':draw_edge, 'genome':genome,'nnames':nnames,'hnGroup':assInfo['hnGroup'],'hLinks':assInfo['hLinks']}
        

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
            lenList.append(arr[1])
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
    #species = para.get("species")
    nodeFile = os.path.join(BinDir,"upload","node.info")
    edgeFile = os.path.join(BinDir,"upload","edge.info")
    chrListFile = os.path.join(BinDir,"upload","chr.list")
    assListFile = os.path.join(BinDir,"upload","ass.list")
    pathFile = os.path.join(BinDir,"upload","path.info")
    sepFile = os.path.join(BinDir,"upload","sep.info")
    
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
    graphInfo = formatGraph(nodeFile,edgeFile,pathFile,sepFile,ass,sChr,sStart,sEnd,ex,wStart,wWidth,wCut,y)
    graphInfo.update(chrList)
    graphInfo['ass'] = assList
    #
    return JsonResponse(graphInfo)

@csrf_exempt
def initGraph(request):
    para = request.POST
    #species = para.get("species")
    nodeFile = os.path.join(BinDir,"upload","node.info")
    edgeFile = os.path.join(BinDir,"upload","edge.info")
    chrListFile = os.path.join(BinDir,"upload","chr.list")
    assListFile = os.path.join(BinDir,"upload","ass.list")
    pathFile = os.path.join(BinDir,"upload","path.info")
    sepFile = os.path.join(BinDir,"upload","sep.info")
    
    chrList = readChr(chrListFile)
    assList = readAss(assListFile)
    
    sChr = chrList["nameList"][0]
    sStart = 1
    sEnd = 10000
    
    ex = 1000000
    wStart = 50
    wWidth = 800
    wCut = 2000
    y = 300
    
    ass = "0"
    graphInfo = formatGraph(nodeFile,edgeFile,pathFile,sepFile,ass,sChr,sStart,sEnd,ex,wStart,wWidth,wCut,y)
    graphInfo.update(chrList)
    graphInfo['ass'] = assList
    
    return JsonResponse(graphInfo)    
    
@csrf_exempt            
def searchNode(request):
    para = request.POST
    #species = para.get("species")
    node = para.get('seg')
    covFile = os.path.join(BinDir,"upload","cover.info")
    nodeFile = os.path.join(BinDir,"upload","node.info")
    sepFile = os.path.join(BinDir,"upload","sep.info")
    sep = getSep(sepFile)
    header = []
    cov = []
    if os.path.exists(covFile):
        with open(covFile) as cf:
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
                nodeAss = "REF" if len(nodeAssArr) == 1 else nodeAssArr[0] + sep + nodeAssArr[1]
                nodeChr = nodeAssArr[-1]
                nodeStart = arr[2]
                nodeEnd = arr[3]
                break
    return JsonResponse({'ass':header,'cov':cov,'nodeAss':nodeAss,'nodeChr':nodeChr,'nodeStart':nodeStart,'nodeEnd':nodeEnd})   

 




# *****************************************************************************

import io

import numpy as np

# *****************************************************************************


def s(f):

    return ('"{:05.1f}"'.format(f))

# *****************************************************************************


def sn(f):

    return ('{:05.1f}'.format(f))

# *****************************************************************************


def svg_defaultgrid(inf_dat, f_lab):

    infDat = np.array(inf_dat)
    ngrid03 = int(np.sqrt(infDat.shape[0]))
    infect = (np.sum(infDat, axis=0)).tolist()

    # Color definitions
    blk = '"rgb(  0,  0,  0)"'
    gry = '"rgb(220,220,220)"'
    # blu = '"rgb(  0,114,189)"'
    red = '"rgb(217, 83, 25)"'
    grn = '"rgb(  0,127,  0)"'
    wht = '"rgb(255,255,255)"'

    # Sizes
    gxtot = 300.0
    gytot = 300.0
    margin = 30.0
    gsep = 40.0
    ylabf = 20.0
    yunitf = 16.0
    xlabf = 20.0
    xunitf = 16.0
    ctoty = 50.0
    posmark = 4.0
    fRate = 50.0

    maxinf = int(np.ceil(1.2*np.max(infect)))
    maxinf = maxinf+4 - maxinf % 4
    infect = [round(float(val)/float(maxinf), 3) for val in infect]

    tmp_yvec = [int(4*maxinf/4), int(3*maxinf/4),
                int(2*maxinf/4), int(1*maxinf/4)]
    tmp_ylab = list()
    for kv in tmp_yvec:
        if (kv < 1000):
            tmp_ylab.append('{:d}'.format(kv))
        elif (kv < 10000):
            tmp_ylab.append('{:3.1f}k'.format(kv/1000))
        else:
            tmp_ylab.append('{:d}k'.format(int(kv/1000)))

    # Labels
    # Automatically overwritten; must be odd
    g1xLabL = ['XX', 'XX', 'XX', 'XX', 'XX', 'XX', 'XX']
    g1yLabL = [tmp_ylab[0], '', tmp_ylab[1], '',
               tmp_ylab[2], '', tmp_ylab[3], '', '0']

    g1xLabNum = len(g1xLabL)
    g1yLabNum = len(g1yLabL)

    gridmult = 10
    tmodval = 10
    nibars = (g1xLabNum-1)*tmodval
    t_start = 0
    t_end = t_start+int(np.ceil(gridmult*len(infect)/tmodval))

    # Prep calcs
    ylab = 1.5*ylabf
    yunit = 1.5*yunitf
    xlab = 1.5*xlabf
    xunit = 1.5*xunitf

    xtot = 2.0*margin + 2.0*gxtot + gsep + 2.0*ylab + 2.0*yunit
    ytot = 2.0*margin + 1.0*gytot + 0.00 + 1.0*xlab + 1.0*xunit
    xall = xtot
    yall = ytot + ctoty

    g1xLaby = margin+gytot+xunit+xlab
    g1xUnity = margin+gytot+xunit
    g1yLabx = margin+ylab/3.0
    g1yUnitx = margin+ylab+yunitf/3.0

    # g3xLaby = margin+gytot+xunit+xlab
    # g3xUnity = margin+gytot+xunit
    # g3yLabx = margin+ylab+yunit+gxtot+gsep+ylab/2.0
    # g3yUnitx = margin+ylab+yunit+gxtot+gsep+ylab+yunitf/3.0

    g1minx = margin+ylab+yunit
    g1maxx = margin+ylab+yunit+gxtot
    g1totx = gxtot
    g1miny = margin
    g1maxy = margin+gytot
    g1toty = gytot

    g3minx = margin+ylab+yunit+gxtot+gsep+ylab+yunit
    # g3maxx = margin+ylab+yunit+gxtot+gsep+ylab+yunit+gxtot
    # g3totx = gxtot
    g3miny = margin
    # g3maxy = margin+gytot
    # g3toty = gytot

    sq_scale03 = gytot/ngrid03
    xpad03 = (gxtot-gytot)/2.0

    cminy = ytot

    c1totx = ctoty
    c2totx = 0.90*gxtot
    c3totx = ctoty

    cminx = xall/2.0 - (c1totx+c2totx+c3totx)/2.0
    c1minx = cminx
    c2minx = c1minx + c1totx
    c3minx = c2minx + c2totx

    fid01 = io.StringIO()

    # Header info
    fid01.write('<?xml version="1.0" encoding="UTF-8"?>'+'\n')
    fid01.write(''+'\n')
    fid01.write('<svg xmlns="http://www.w3.org/2000/svg"'+'\n')
    fid01.write('     width='+s(xall)+'\n')
    fid01.write('     height='+s(yall)+'\n')
    fid01.write('     viewbox="0 0 '+sn(xall)+' '+sn(yall)+'"'+'\n')
    fid01.write('     id="SVGDoc"\n')
    fid01.write('     onload="Init()">'+'\n')
    fid01.write(''+'\n')
    fid01.write('  <title>IDM Simulation Example</title>'+'\n')
    fid01.write(''+'\n')
    fid01.write('  <desc></desc>'+'\n')
    fid01.write(''+'\n')
    fid01.write('  <script language="JavaScript">'+'\n')
    fid01.write('    <![CDATA['+'\n')
    fid01.write('      var SVGDoc = null;'+'\n')
    fid01.write('      var intervalHandle = null;'+'\n')
    fid01.write('      var infDat = '+str(infect)+';'+'\n')
    fid01.write('      var infBlk = '+str(infDat.tolist())+';'+'\n')
    fid01.write('      var tVal = 0;'+'\n')
    fid01.write('      var junk = 0;'+'\n')
    fid01.write('      var gInfMax = '+str(np.amax(infDat))+';'+'\n')
    fid01.write('      var loopIt = 0;'+'\n')
    fid01.write('      var playIt = 0;'+'\n')
    fid01.write('      var pauseIt = 0;'+'\n')
    fid01.write('      var yInfB = '+sn(g1maxy)+';'+'\n')
    fid01.write('      var yInfT = '+sn(g1toty)+';'+'\n')
    fid01.write('      var year0 = '+str(t_start)+';'+'\n')
    fid01.write('      var yearL = '+str(t_end)+';'+'\n')
    fid01.write('      var posmin = '+str(c2minx+0.05*c2totx)+';'+'\n')
    fid01.write('      var posmax = '+str(c2minx+0.95*c2totx)+';'+'\n')
    fid01.write('      var posmark = '+str(posmark)+';'+'\n')
    fid01.write(''+'\n')
    fid01.write('      function Init()'+'\n')
    fid01.write('      {'+'\n')
    fid01.write('        SVGDoc = document.getElementById("SVGDoc");'+'\n')
    fid01.write('        rePaint();'+'\n')
    fid01.write('      };'+'\n')
    fid01.write(''+'\n')
    fid01.write('      function ctrlOn(event)'+'\n')
    fid01.write('      {'+'\n')
    fid01.write('        var xVal = event.clientX;'+'\n')
    fid01.write('        if(playIt == 1)'+'\n')
    fid01.write('        {'+'\n')
    fid01.write('          pauseIt = 1;'+'\n')
    fid01.write('          clearInterval(intervalHandle);'+'\n')
    fid01.write('        }'+'\n')
    fid01.write('        movPosbar(xVal);'+'\n')
    fid01.write('        SVGDoc.addEventListener("mousemove",dragPosbar);'+'\n')
    fid01.write('      };'+'\n')
    fid01.write(''+'\n')
    fid01.write('      function ctrlOff()'+'\n')
    fid01.write('      {'+'\n')
    fid01.write('        if(playIt == 1 && pauseIt ==1)'+'\n')
    fid01.write('        {'+'\n')
    fid01.write('          pauseIt = 0;'+'\n')
    fid01.write('          intervalHandle = setInterval(goTime,'+str(fRate)+');'+'\n')
    fid01.write('        }'+'\n')
    fid01.write('        SVGDoc.removeEventListener("mousemove",dragPosbar);'+'\n')
    fid01.write('      };'+'\n')
    fid01.write(''+'\n')
    fid01.write('      function goTime()'+'\n')
    fid01.write('      {'+'\n')
    fid01.write('        tVal++;'+'\n')
    fid01.write('        if(tVal < infDat.length && tVal >= 0)'+'\n')
    fid01.write('        {'+'\n')
    fid01.write('          rePaint();'+'\n')
    fid01.write('        }'+'\n')
    fid01.write('        if(loopIt == 1 && tVal >= (infDat.length+10))'+'\n')
    fid01.write('        {'+'\n')
    fid01.write('          tVal = -10;'+'\n')
    fid01.write('        }'+'\n')
    fid01.write('        if(tVal >= (infDat.length+20))'+'\n')
    fid01.write('        {'+'\n')
    fid01.write('          tVal = (infDat.length+20);'+'\n')
    fid01.write('        }'+'\n')
    fid01.write('      };'+'\n')
    fid01.write(''+'\n')
    fid01.write('      function togRep()'+'\n')
    fid01.write('      {'+'\n')
    fid01.write('        if(loopIt == 0)'+'\n')
    fid01.write('        {'+'\n')
    fid01.write('          loopIt = 1;'+'\n')
    fid01.write('          var bID = SVGDoc.getElementById("repeatshape");'+'\n')
    fid01.write('          bID.setAttribute("stroke",'+grn+');'+'\n')
    fid01.write('        }'+'\n')
    fid01.write('        else'+'\n')
    fid01.write('        {'+'\n')
    fid01.write('          loopIt = 0;'+'\n')
    fid01.write('          var bID = SVGDoc.getElementById("repeatshape");'+'\n')
    fid01.write('          bID.setAttribute("stroke",'+blk+');'+'\n')
    fid01.write('        }'+'\n')
    fid01.write('      };'+'\n')
    fid01.write(''+'\n')
    fid01.write('      function togPlay()'+'\n')
    fid01.write('      {'+'\n')
    fid01.write('        if(playIt == 0)'+'\n')
    fid01.write('        {'+'\n')
    fid01.write('          playIt = 1;'+'\n')
    fid01.write('          var bID = SVGDoc.getElementById("playshape");'+'\n')
    fid01.write('          bID.setAttribute("fill",'+grn+');'+'\n')
    fid01.write('          intervalHandle = setInterval(goTime,'+str(fRate)+');'+'\n')
    fid01.write('        }'+'\n')
    fid01.write('        else'+'\n')
    fid01.write('        {'+'\n')
    fid01.write('          playIt = 0;'+'\n')
    fid01.write('          var bID = SVGDoc.getElementById("playshape");'+'\n')
    fid01.write('          bID.setAttribute("fill",'+wht+');'+'\n')
    fid01.write('          clearInterval(intervalHandle);'+'\n')
    fid01.write('        }'+'\n')
    fid01.write('      };'+'\n')
    fid01.write(''+'\n')
    fid01.write('      function rePaint()'+'\n')
    fid01.write('      {'+'\n')
    fid01.write('        paintInfBar();'+'\n')
    fid01.write('        paintDemoGrid();'+'\n')
    fid01.write('        paintInfAxis();'+'\n')
    fid01.write('        paintPosbar();'+'\n')
    fid01.write('      };'+'\n')
    fid01.write(''+'\n')
    fid01.write('      function paintDemoGrid()'+'\n')
    fid01.write('      {'+'\n')
    fid01.write('        for(k1 = 0; k1 < infBlk.length; k1++)'+'\n')
    fid01.write('        {'+'\n')
    fid01.write('          var bID = SVGDoc.getElementById("gridR"+(k1).toString())'+'\n')
    fid01.write('          bID.setAttribute("fill-opacity",(infBlk[k1][tVal]/gInfMax).toString());'+'\n')
    fid01.write('        }'+'\n')
    fid01.write('      };'+'\n')
    fid01.write(''+'\n')
    fid01.write('      function paintInfBar()'+'\n')
    fid01.write('      {'+'\n')
    fid01.write('        for(k1 = 0; k1 < '+str(nibars)+'; k1++)'+'\n')
    fid01.write('        {'+'\n')
    fid01.write('          var bID = SVGDoc.getElementById("infR"+(k1).toString())'+'\n')
    fid01.write('          var dVal = tVal+k1-'+str(nibars/2)+';'+'\n')
    fid01.write('          if(dVal >= 0 && dVal < infDat.length)'+'\n')
    fid01.write('          {'+'\n')
    fid01.write('            bID.setAttribute("y",(yInfB-yInfT*infDat[dVal]).toString());'+'\n')
    fid01.write('            bID.setAttribute("height",(yInfT*infDat[dVal]).toString());'+'\n')
    fid01.write('          }'+'\n')
    fid01.write('          else'+'\n')
    fid01.write('          {'+'\n')
    fid01.write('            bID.setAttribute("y",(0).toString());'+'\n')
    fid01.write('            bID.setAttribute("height",(0).toString());'+'\n')
    fid01.write('          }'+'\n')
    fid01.write('        }'+'\n')
    fid01.write('      };'+'\n')
    fid01.write(''+'\n')
    fid01.write('      function paintInfAxis()'+'\n')
    fid01.write('      {'+'\n')
    fid01.write('        for(k1 = 0; k1 < '+str(g1xLabNum)+'; k1++)'+'\n')
    fid01.write('        {'+'\n')
    fid01.write('          var bID1 = SVGDoc.getElementById("yrlab"+(k1).toString());'+'\n')
    fid01.write('          var bID2 = SVGDoc.getElementById("yrgrid"+(k1).toString());'+'\n')
    fid01.write('          var tMod = tVal%'+str(tmodval)+';'+'\n')
    fid01.write('          var dVal = '+str(gridmult)+'*((tVal-tMod)/'+str(tmodval)+'-'+str((g1xLabNum-1)/2)+'+k1)+year0;'+'\n')
    fid01.write('          var newX = '+str(g1minx)+'+(k1-tMod/'+str(tmodval)+')*'+str(g1totx)+'/'+str(g1xLabNum-1)+';'+'\n')
    fid01.write('          if(dVal >= year0 && dVal <= yearL)'+'\n')
    fid01.write('          {'+'\n')
    fid01.write('            bID1.textContent = dVal.toString();'+'\n')
    fid01.write('          }'+'\n')
    fid01.write('          else'+'\n')
    fid01.write('          {'+'\n')
    fid01.write('            bID1.textContent = "";'+'\n')
    fid01.write('          }'+'\n')
    fid01.write('          if(newX >= '+str(g1minx)+')'+'\n')
    fid01.write('          {'+'\n')
    fid01.write('            bID1.setAttribute("x",newX.toString());'+'\n')
    fid01.write('            bID2.setAttribute("x1",newX.toString());'+'\n')
    fid01.write('            bID2.setAttribute("x2",newX.toString());'+'\n')
    fid01.write('          }'+'\n')
    fid01.write('          else'+'\n')
    fid01.write('          {'+'\n')
    fid01.write('            bID1.textContent = "";'+'\n')
    fid01.write('            bID2.setAttribute("x1",'+str(g1minx)+');'+'\n')
    fid01.write('            bID2.setAttribute("x2",'+str(g1minx)+');'+'\n')
    fid01.write('          }'+'\n')
    fid01.write('        }'+'\n')
    fid01.write('      };'+'\n')
    fid01.write(''+'\n')
    fid01.write('      function paintPosbar()'+'\n')
    fid01.write('      {'+'\n')
    fid01.write('        var bID = SVGDoc.getElementById("posMark");'+'\n')
    fid01.write('        var fracpos;'+'\n')
    fid01.write('        if(tVal < 0)'+'\n')
    fid01.write('        {'+'\n')
    fid01.write('          fracpos = 0.0;'+'\n')
    fid01.write('        }'+'\n')
    fid01.write('        else if(tVal > (infDat.length-1))'+'\n')
    fid01.write('        {'+'\n')
    fid01.write('          fracpos = 1.0;'+'\n')
    fid01.write('        }'+'\n')
    fid01.write('        else if(infDat.length == 1)'+'\n')
    fid01.write('        {'+'\n')
    fid01.write('          fracpos = 0.5;'+'\n')
    fid01.write('        }'+'\n')
    fid01.write('        else'+'\n')
    fid01.write('        {'+'\n')
    fid01.write('          fracpos = tVal/(infDat.length-1);'+'\n')
    fid01.write('        }'+'\n')
    fid01.write('        bID.setAttribute("x",posmin+(posmax-posmin)*fracpos);'+'\n')
    fid01.write('      };'+'\n')
    fid01.write(''+'\n')
    fid01.write('      function dragPosbar(event)'+'\n')
    fid01.write('      {'+'\n')
    fid01.write('        var xVal = event.clientX;'+'\n')
    fid01.write('        movPosbar(xVal);'+'\n')
    fid01.write('      };'+'\n')
    fid01.write(''+'\n')
    fid01.write('      function movPosbar(xVal)'+'\n')
    fid01.write('      {'+'\n')
    fid01.write('        nTime = Math.round(infDat.length*(xVal-posmin)/(posmax-posmin));'+'\n')
    fid01.write('        if(nTime < 0)'+'\n')
    fid01.write('        {'+'\n')
    fid01.write('          nTime = 0;'+'\n')
    fid01.write('        }'+'\n')
    fid01.write('        else if(nTime > (infDat.length-1))'+'\n')
    fid01.write('        {'+'\n')
    fid01.write('          nTime = (infDat.length-1);'+'\n')
    fid01.write('        }'+'\n')
    fid01.write('        tVal = nTime;'+'\n')
    fid01.write('        rePaint();'+'\n')
    fid01.write('      };'+'\n')
    fid01.write(''+'\n')
    fid01.write('    ]]>'+'\n')
    fid01.write('  </script>'+'\n')

    # Controls
    fid01.write('  <!--  Controls Box -->'+'\n')
    fid01.write('  <g stroke='+blk+' stroke-width="3.0">'+'\n')
    fid01.write('    <path d="M'+sn(c1minx+0.25*c1totx)+' '+sn(cminy+0.2*ctoty)+' L'+sn(c1minx+0.25*c1totx)+' '+sn(cminy+0.8*ctoty)+' L'+sn(c1minx+0.8*c1totx)+' '+sn(cminy+0.5*ctoty)+' Z" fill='+wht+' id="playshape"/>'+'\n')
    fid01.write('    <path d="M'+sn(c3minx+0.2*c3totx)+' '+sn(cminy+0.5*ctoty)+' A'+sn(0.3*c3totx)+' '+sn(0.3*ctoty)+' 0.0 1 0 '+sn(c3minx+0.5*c3totx)+' '+sn(cminy+0.2*ctoty)+' M'+sn(c3minx+0.5*c3totx)+' '+sn(cminy+0.2*ctoty)+' L'+sn(c3minx+0.60*c3totx)+' '+sn(cminy+0.3*ctoty)+' M'+sn(c3minx+0.5*c3totx)+' '+sn(cminy+0.2*ctoty)+' L'+sn(c3minx+0.6*c3totx)+' '+sn(cminy+0.1*ctoty)+'" fill='+wht+' stroke-linecap="round" id="repeatshape"/>'+'\n')
    fid01.write('    <rect x='+s(c2minx+0.05*c2totx)+' y='+s(cminy+0.45*ctoty)+' width='+s(0.9*c2totx)+' height='+s(0.10*ctoty)+' stroke='+blk+' fill="transparent" />'+'\n')
    fid01.write('    <rect x='+s(c2minx+0.05*c2totx-posmark/2.0)+' y='+s(cminy+0.35*ctoty)+' width='+s(posmark)+' height='+s(0.30*ctoty)+' stroke='+blk+' fill='+blk+' id="posMark"/>'+'\n')
    fid01.write('    <rect x='+s(c2minx)+' y='+s(cminy)+' width='+s(c2totx)+' height='+s(ctoty)+' stroke-width="0.0" fill="transparent" onmousedown="ctrlOn(evt)" onmouseup="ctrlOff()" onmouseout="ctrlOff()"/>'+'\n')
    fid01.write('    <circle cx='+s(c1minx+0.5*c1totx)+' cy='+s(cminy+0.5*ctoty)+' r='+s(0.5*ctoty)+' stroke-width="0.0" fill="transparent" onclick="togPlay()"/>'+'\n')
    fid01.write('    <circle cx='+s(c3minx+0.5*c3totx)+' cy='+s(cminy+0.5*ctoty)+' r='+s(0.5*ctoty)+' stroke-width="0.0" fill="transparent" onclick="togRep()"/>'+'\n')
    fid01.write('  </g>'+'\n')
    fid01.write(''+'\n')

    # Graph 1 grid
    fid01.write('  <!--  Graph 1 Grid -->'+'\n')
    fid01.write('  <g stroke='+gry+' stroke-width="1.0">'+'\n')
    for k1 in range(g1yLabNum):
        fid01.write('    <line x1='+s(g1minx)+' y1='+s(g1miny+float(k1)*g1toty/float(g1yLabNum-1))+' x2='+s(g1maxx)+' y2='+s(g1miny+float(k1)*g1toty/float(g1yLabNum-1))+'/>'+'\n')

    fid01.write(''+'\n')
    for k1 in range(g1xLabNum):
        fid01.write('    <line x1='+s(g1minx+float(k1)*g1totx/float(g1xLabNum-1))+' y1='+s(g1miny)+' x2='+s(g1minx+float(k1)*g1totx/float(g1xLabNum-1))+' y2='+s(g1maxy)+' id="yrgrid'+str(k1)+'"/>'+'\n')

    fid01.write('  </g>'+'\n')
    fid01.write(''+'\n')

    # Graph 1 bars
    fid01.write('  <!--  Graph 1 Bars -->'+'\n')
    fid01.write('  <g fill='+red+' stroke-width="0.0" stroke-opacity="0.0">'+'\n')
    for k1 in range(nibars):
        fid01.write('    <rect x='+s(g1minx+float(k1)*g1totx/float(nibars))+' y='+s(g1maxy)+' width='+s(g1totx/float(nibars))+' height='+s(0)+' id="infR{:d}"'.format(k1)+'/>'+'\n')

    fid01.write('  </g>'+'\n')
    fid01.write(''+'\n')

    # Grid 3 patches
    fid01.write('  <!--  Grid 3 Patches -->'+'\n')
    fid01.write('  <g stroke='+blk+' fill='+red+' fill-opacity="0.0" stroke-width="0.3" stroke-opacity="1.0">'+'\n')
    for k1 in range(ngrid03*ngrid03):
        xv = g3minx + xpad03 + sq_scale03*(k1 % ngrid03)
        yv = g3miny + sq_scale03*(ngrid03 - k1//ngrid03 - 1)
        fid01.write('    <rect x='+s(xv)+' y='+s(yv)+' width='+s(sq_scale03)+' height='+s(sq_scale03)+' id="gridR{:d}"'.format(k1)+'/>'+'\n')

    fid01.write('  </g>'+'\n')
    fid01.write(''+'\n')

    # Graph 1 divider
    fid01.write('  <!--  Graph 1 Divider -->'+'\n')
    fid01.write('  <g stroke='+blk+' stroke-width="1.5" stroke-dasharray="5,5">'+'\n')
    fid01.write('    <line x1='+s(g1minx+g1totx/2.0)+' y1='+s(g1miny)+' x2='+s(g1minx+g1totx/2.0)+' y2='+s(g1maxy)+'/>'+'\n')
    fid01.write('  </g>'+'\n')
    fid01.write(''+'\n')

    # Graph 1 border
    fid01.write('  <!--  Graph 1 Border -->'+'\n')
    fid01.write('  <g stroke='+blk+' stroke-width="2.0">'+'\n')
    fid01.write('    <line x1='+s(g1minx)+' y1='+s(g1miny)+' x2='+s(g1maxx)+' y2='+s(g1miny)+'/>'+'\n')
    fid01.write('    <line x1='+s(g1minx)+' y1='+s(g1maxy)+' x2='+s(g1maxx)+' y2='+s(g1maxy)+'/>'+'\n')
    fid01.write('    <line x1='+s(g1minx)+' y1='+s(g1miny)+' x2='+s(g1minx)+' y2='+s(g1maxy)+'/>'+'\n')
    fid01.write('    <line x1='+s(g1maxx)+' y1='+s(g1miny)+' x2='+s(g1maxx)+' y2='+s(g1maxy)+'/>'+'\n')
    fid01.write('  </g>'+'\n')
    fid01.write(''+'\n')

    # Graph 1 labels
    fid01.write('  <!--  Graph 1 Labels -->'+'\n')
    fid01.write('  <g stroke='+blk+' text-anchor="middle" font-family="Sans,Arial">'+'\n')
    fid01.write('    <text x='+s(g1minx+g1totx/2.0)+' y='+s(g1xLaby)+' font-size='+s(xlabf)+'>Time</text>'+'\n')
    for k1 in range(g1xLabNum):
        fid01.write('    <text x='+s(g1minx+float(k1)*g1totx/float(g1xLabNum-1))+' y='+s(g1xUnity)+' font-size='+s(xunitf)+' id="yrlab'+str(k1)+'">'+g1xLabL[k1]+'</text>'+'\n')

    fid01.write('    <text x='+s(g1yLabx)+' y='+s(g1miny+g1toty/2.0)+' font-size='+s(ylabf)+' transform="rotate(-90.0 '+sn(g1yLabx)+','+sn(g1miny+g1toty/2.0)+')">Incidence</text>'+'\n')
    for k1 in range(g1yLabNum):
        fid01.write('    <text x='+s(g1yUnitx)+' y='+s(g1miny+float(k1)*g1toty/float(g1yLabNum-1)+yunitf/3.0)+' font-size='+s(yunitf)+'>'+g1yLabL[k1]+'</text>'+'\n')

    fid01.write('  </g>'+'\n')
    fid01.write(''+'\n')

    # End SVG
    fid01.write(''+'\n')
    fid01.write('</svg>'+'\n')

    # Write file
    with open('networkviz_01_{:s}.svg'.format(f_lab), 'w') as dump01:
        dump01.write(fid01.getvalue())

    return None

# *****************************************************************************

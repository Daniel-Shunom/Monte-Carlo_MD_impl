"""
.. module::  MCTree.py
   :platform: Unix, Windows, Linux
   :synopsis: Primary code responsible for upper level mangagement of the tree data structure.

.. moduleauthor:: Troy D Loeffler <tloeffler@anl.gov>


"""

from random import choice, random, shuffle
from math import sqrt, log
from Node import Node
#=================================================
class InitialDataNotGiven(Exception):
    '''
     The user forgot to pass in a Data Object to kick things off.
    '''
    pass
#=================================================
class UnnormalizedProb(Exception):
    pass
#=================================================
class MissingHeadNode(Exception):
    '''
     HeadNode is somehow not present in the restart file
    '''
    pass
#=================================================
class NodeIDNotFound(Exception):
    '''
     Somehow in the restart file a node was referenced in playout section
     and yet was not defined in the main section
    '''
    pass

#=================================================
class BadFormating(Exception):
    '''
     Restart file format does not match what is expected
    '''
    pass

#=================================================
class UnknownFormat(Exception):
    '''
     Restart file format does not match what is expected
    '''
    pass
#=================================================
def randomlist(inlist, problist):
    sumint = 0.0
    indx = -1
    randNum = random()
    while sumint < randNum:
        indx += 1
        sumint += problist[indx] 
#        print sumint, randNum, indx

    return inlist[indx]

#=================================================
def DefaultPath(nodelist):
    print("Selecting Pathway")
    def PathScore(node):
        score = node.getbestplayscore()
        score2 = node.getbestscore()
        scoremin = min(score, score2)
        print(node.getid(), score, score2, scoremin)
        return scoremin
    selection = sorted(nodelist, key=lambda x:PathScore(x))[0]
    print("Node %s is now the new headnode"%(selection.getid()))
    return selection, -1
#=================================================
def DefaultPrune(nodelist, maxdepth=3):
    print("Pruning Tree")
    def PathScore(node):
        score = node.getbestplayscore()
        print(node.getid(), score)
        return score
    shortlist = [x for x in nodelist  if x.getdepth() < maxdepth ]
    selection = sorted(shortlist, key=lambda x:PathScore(x))[-1]
    print("Node %s and it's child nodes have been pruned"%(selection.getid()))
    return [selection]

#=================================================
class Tree(object):
    '''
     Primary Tree Object which manages all opperations that will be performed. 
    '''
    #--------------------------------------------
    def __init__(self, seeddata=None, selectfunction=None, prunefunction=None,pathfunction=None, playchoice=None, playouts=2, headexpansion=5):
        if seeddata is None:
            raise InitialDataNotGiven("The seeddata argument must be specified")
        self.headnode = Node(data=seeddata)
        self.headnode.setid(0)
        self.headexpansion = headexpansion
        self.curid = 0
        self.playouts = playouts
        self.exploreconst = sqrt(2.0)
        if selectfunction is None:
            self.select = UCB
        else:
            self.select = selectfunction

        if prunefunction is None:
            self.prune = DefaultPrune
        else:
            self.prune = prunefunction

        if pathfunction is None:
            self.pathfunction = DefaultPath
        else:
            self.pathfunction = pathfunction

        self.playselect = None

    #--------------------------------------------
    def __str__(self):
        minscore = min(self.headnode.getbestscore(), self.headnode.getbestplayscore())
        outstr = "Best Score: %s"%(minscore)
        return outstr
    #--------------------------------------------
    def expand(self, nExpansions=5, writeevery=None, saveevery=None):
        '''
         This function adds a child node to the selected node in the tree.  This is the primary
         funcion for increasing the tree size. 


        '''

        if self.headnode.getvisits() == 0:
            if self.playouts > 0:
                self.headnode.simulate(self.playouts)
#                self.headnode.backpropigate()


        while len(self.headnode.getchildren()) < self.headexpansion:
            self.curid += 1
            newnode = self.headnode.spawnchild(self.curid)
            if newnode is None:
                print("Node is degenerate, trying another...")
                continue

            newnode.evaluate()
#            newnode.associateplayouts()
            if self.playouts > 0:
                newnode.simulate(self.playouts)
#            newnode.backpropigate()

        self.writenewick("treestructure.nh")
        self.writeanalytics("analytics.dat")
        for isearch in range(nExpansions):
            nodelist = self.headnode.getnodelist()
            targnode = self.select(nodelist, self.exploreconst)
            self.curid += 1
            newnode = targnode.spawnchild(self.curid)
            if newnode is None:
                print("Node is degenerate, trying another...")
                continue
            newnode.evaluate()
#            newnode.associateplayouts()
            if self.playouts > 0:
                newnode.simulate(self.playouts)
#                newnode.backpropigate()

            if writeevery is not None:
                if isearch%writeevery == 0:
                    self.writenewick("treestructure.nh")
                    self.writeanalytics("analytics.dat")

            if saveevery is not None:
                if isearch%saveevery == 0:
                    self.savetree("tree.chkpnt")
    #--------------------------------------------
    def playexpand(self, nExpansions=5, recalc=False):
        '''
         .. note:
         Unlike regular expansion where a brandnew point is created this function
         uses one of the playouts of a selected node to create a new child node.

         Note: Recalc must be set to True if the playout objective and true objective functions
               are different.  This forces the tree to recompute the score using the exact objective function.
        '''

        for isearch in range(nExpansions):
            nodelist = self.headnode.getnodelist()
            targnode = self.select(nodelist, self.exploreconst)
            self.curid += 1
            newnode = targnode.spawnfromplayout(self.curid, self.playselect)
            if newnode is None:
                print("Node is degenerate, trying another...")
                return
            if recalc:
                newnode.evaluate()
#            newnode.associateplayouts()
            if self.playouts > 0:
                newnode.simulate(self.playouts)
    #--------------------------------------------
    def minexpand(self, nExpansions=1, recalc=False, scaleconst=None, mindepth=0):
        '''
          Creates a new node by running local optimization on the parent node and uses the final solution
          as the structure for the new node.

        '''

        if scaleconst is not None:
            self.scaleconstant(scaleconst)
        for isearch in range(nExpansions):
            nodelist = self.headnode.getnodelist()
            nodelist = [x for x in nodelist if x.getdepth() > mindepth]
            if len(nodelist) < 1:
                print("Nothing to Minimize Yet. Skipping for now.")
                return
            nodelist = [x for x in nodelist if not x.isminned()]
            targnode = self.select(nodelist, self.exploreconst)
            if targnode.isminned():
                continue
            self.curid += 1
            newnode, score = targnode.spawnfromminimize(self.curid)
            if newnode is None:
                print("Node is degenerate, trying another...")
                return
            newnode.setscore(score)
#            newnode.associateplayouts()
            if self.playouts > 0:
                newnode.simulate(self.playouts)

        if scaleconst is not None:
            self.scaleconstant(1.0/scaleconst)
    #--------------------------------------------
    def crossover(self, nExpansions=5):
        '''
         Crossover opperation from Genetic Algorithms. This creates
         a new node by selecting two parent nodes to make an offspring state.

        '''
        for isearch in range(nExpansions):
            nodelist = self.headnode.getnodelist()
            targnode = self.select(nodelist, self.exploreconst)
            nodelist = [x for x in self.headnode.getnodelist() if x.getid() != targnode.getid()]
            cnt = 0
#            partnernode = self.select(nodelist, self.exploreconst)
            partnernode = self.select(nodelist, self.exploreconst*1e-1)

            self.curid += 1
            newnode = targnode.spawnfromcrossover(partnernode=partnernode, curid=self.curid)
            if newnode is None:
                print("Node is degenerate, trying another...")
                continue
            newnode.evaluate()
#            newnode.associateplayouts()
            if self.playouts > 0:
                newnode.simulate(self.playouts)

    #--------------------------------------------
    def graft(self, nGraft=1):
        '''
        Takes a node located lower in the tree, disconnects it from it's parent, and
        then ties it to the head node.  This is useful for mixing up the tree structure when
        a window scaling scheme is used.  A node which is grafted will take it's children
        nodes with it and be attached to the head node.

        '''
        for isearch in range(nGraft):
            nodelist = [x for x in self.headnode.getnodelist() if x.getdepth() > 1]
            if len(nodelist) < 1:
                return
            sortlist = sorted(nodelist, key=lambda x: (-x.getbestscore(), x.getdepth()) )
            selection = sortlist[-1]

            selparent = selection.getparent()
            selparent.dissociatechild(selection)
            self.headnode.addchild(selection)
            self.headnode.createdepth()
            print("Node %s has been grafted from Node %s to Node %s"%(selection.getid(), selparent.getid(), self.headnode.getid()))
#            sleep(5)




    #--------------------------------------------
    def simulate(self, nSimulations=5):
        '''
         This is the roll out function. The nSimulation parameter will determine how many nodes 
         are selected for a playout during this function call. 
        '''
        if self.playouts < 1:
            return
        for iSearch in range(nSimulations):
            nodelist = self.headnode.getnodelist()
            targnode = self.select(nodelist, self.exploreconst)
            targnode.simulate(self.playouts)
#            targnode.backpropigate()
    #--------------------------------------------
    def selectpath(self):
        '''
           Advances the tree by selecting a child node to become the new headnode.
           The selection function used requires that the new headnode be returned
           in addition to the delta which is used to readjust the depth variable
           for all the nodes on this given branch. 
        '''

        childlist = self.headnode.getchildren()
        self.headnode, delta = self.pathfunction(childlist)
        self.headnode.setparent(None)
        self.headnode.changedepth(delta=delta)

    #--------------------------------------------
    def prunetree(self, nPrune=1):
        '''
           Used to trim nodes off the tree that are not useful. This is done to clean up memory
           as the tree grows large.
        '''

        for iPrune in range(nPrune):
            nodelist = self.headnode.getchildren()
            prunelist = self.prune(nodelist)
            for node in prunelist:
                if node.getstatus() != "head":
                    parent = node.getparent()
                    parent.dissociatechild(node)

    #--------------------------------------------
    def writenewick(self, outname):
        '''
         Writes the newick tree representation to a file.  This can be used
         by visualization software such as Ete3 to create human readable tree diagrams
        '''
        print("Writing Tree Analyitics")
        newwick = self.headnode.getnewick()
        with open(outname, "a") as outfile:
            outfile.write(newwick+"\n")
    #--------------------------------------------
    def writeanalytics(self, outname):
        '''
         Writes the node analytics to a file. 
        '''
        minscore = min(self.headnode.getbestscore(), self.headnode.getbestplayscore())
        print("Best Score: %s"%(minscore))
        infolist = self.headnode.getratiolist()
        infolist = sorted(infolist, key=lambda x:x[0])
#        formatstr = ' '.join(["%s" for x in ratiolist])
#        formatstr = formatstr + "\n"
        with open(outname, "a") as outfile:
            for line in infolist:
                outfile.write(' '.join([str(x) for x in line]) + "\n")
            outfile.write("\n")

    #--------------------------------------------
    def getdepthhist(self, normalize=False):
        '''
         Returns a histogram showing the distribution of the nodes by tree depth. 
        '''
        nodelist = self.headnode.getnodelist()
        maxdepth = max([x.getdepth() for x in nodelist])
        depthhist = [0 for x in range(maxdepth+1)]

        for node in nodelist:
            depth = node.getdepth()
            depthhist[depth] += 1
        if normalize:
            norm = sum(depthhist)
            depthhist = [float(x)/norm for x in depthhist]

        return depthhist
    #--------------------------------------------
    def getbestnode(self, sortkey=None):
        nodelist = self.headnode.getnodelist()
        if sortkey is None:
            nodelist = sorted(nodelist, key=lambda x:x.getscore())
        else:
            nodelist = sorted(nodelist, key=sortkey)
        bestnode = nodelist[0]
        return bestnode
    #--------------------------------------------
    def getheadexpand(self):
        return self.headexpansion
    #--------------------------------------------
    def setheadexpand(self, newhead):
        self.headexpansion = newhead
        print("Head Expansion is now %s"%(newhead))
    #--------------------------------------------
    def setconstant(self, newconst):
        self.exploreconst = newconst
        print("New Exploration Constant is: %s"%(self.exploreconst))
    #--------------------------------------------
    def scaleconstant(self, scale):
        self.exploreconst *= scale
        print("New Exploration Constant is: %s"%(self.exploreconst))
    #--------------------------------------------
    def getbestscore(self):
        return min(self.headnode.getbestscore(), self.headnode.getbestplayscore())
    #--------------------------------------------
    def getnodelist(self):
        return self.headnode.getnodelist()

    #--------------------------------------------
    def loadtree(self, filename, seeddata):
        '''
         Reloads a tree from a previous run. 
         File Format
           Number of nodes
              NodeID  ParentID  Nodescore | NodeData
              ....
              ....

           Playout start
              Node ID NPlayouts
                PlayScore | PlayData
              Used Playout Numbers
              ....
              ....
        '''
        print("Loading Tree Structure from File: %s"%(filename))
        with open(filename, "r") as loadfile:
            parentdic = {}
            scoredic = {}
            datadic = {}

            try:
                #Collect the number of nodes to make from the first line
                firstline = loadfile.readline()
                nNodes = int(firstline.split()[0])
                headnodeid = None
                for iLine in range(nNodes):
                    line = loadfile.readline()
                    col = line.split("|")
                    col2 = col[0].split()
#                    print(col)
#                    print(col2)
                    try:
                        nodeid, parentid, nodescore = int(col2[0]), int(col2[1]), float(col2[2])
                        nodedata = col[1]
                    except:
                        col2 = col[0].split()
                        nodeid, parentid, nodescore = int(col2[0]), None, float(col2[2])
                        nodedata = col[1]
                        headnodeid = nodeid
                    nodedata.replace("\n","")
                    nodedata = nodedata.strip()
#                    print(nodedata)
                    parentdic[nodeid] = parentid
                    scoredic[nodeid] = nodescore
                    datadic[nodeid] = nodedata

                #Create the node classes that will be used to construct the tree.  They are initially left dissociated.
                tempnodelist = []
                for nodeid in parentdic:
                    nodescore = scoredic[nodeid]
                    nodedata = datadic[nodeid]

                    newdata = seeddata.newdataobject()
                    newnode = Node(data=newdata, computescore=False)
                    print("Node %s Created..."%(nodeid))

                    newnode.setid(nodeid)
                    newnode.setscore(nodescore)
                    newnode.loaddatastr(nodedata)
                    tempnodelist.append(newnode)

                #Begin associating the nodes with their respective parent nodes. This step builds the tree
                #structure.
                headnode = next((x for x in tempnodelist if x.getid() == headnodeid), None)
                if headnode is None:
                    raise MissingHeadNode("The HeadNode was not properly defined in the restart file.")
                self.headnode = headnode
                self.headnode.setparent(None)

                for node in tempnodelist:
                    parentID = parentdic[node.getid()]
                    if parentID is not None:
                        parentnode = next((x for x in tempnodelist if x.getid() == parentID), None)
                        if parentnode is None:
                            raise NodeIDNotFound("A Node's Parent Node was referenced in the restart file, but was not properly defined in the header section of the restart file!")
                        parentnode.addchild(node)

                # Assign proper depth variables for each node.
                nodelist = self.headnode.getnodelist()

                #Begin collecting playout data from the file
                blankCnt = 0
                print("Node Structure Loaded, Loading Playouts.....")
                while True:
                    line = loadfile.readline()
                    if line is None:
                        break #EOF
                    if line == '': 
                        blankCnt += 1
                        if blankCnt > 10:
                            break
                    else:
                        blankCnt = 0
                    if "Node" in line:
                        col = line.split()
                        nodeID, nPlayouts = int(col[1]), int(col[2])
                        currentnode = next((x for x in nodelist if x.getid() == nodeID), None)
                        if currentnode is None:
                            raise NodeIDNotFound("A Node was referenced in the Playout section of the restart file, but it was not found in the main node definition section!")
                        for iLine in range(nPlayouts):
                            readpoint = loadfile.tell()
                            playline = loadfile.readline()
                            if "Node" in playline:
                                raise BadFormating


#                            print(playline)
                            if playline is None:
                                raise
                            playcol = playline.split("|")
                            try:
                                playscore, playdata = float(playcol[0]), playcol[1]

                            except ValueError:
                                print("Error Reading Node %s"%(nodeID))
                                print("Expected Playouts %s"%(nPlayouts))
                                print("Bad Line: %s"%(playline))
                                loadfile.seek(readpoint)
                                break
                            playdata = playdata.strip()
                            playdata = playdata.replace("\n","")
                            currentnode.addplayout(playdata, playscore)
                        usedline = loadfile.readline()
                        usedcol = usedline.split()
                        usedlist = [int(x) for x in usedcol]
                        currentnode.setusedlist(usedlist)

            except MissingHeadNode:
                raise MissingHeadNode("The specified input file did not contain a head node!")

            except UnknownFormat:
                raise UnknownFormat("The file given to loadtree does not match a known MCTS tree format.")
        
        self.headnode.createdepth()
        self.curid = max([x.getid() for x in nodelist]) + 1
        

    #-------------------------------------------------------------------------------------------------------------------
    def savetree(self, filename):
        print("Saving the tree to file: %s"%(filename))
        nodelist = self.headnode.getnodelist()
        nodelist = sorted(nodelist, key= lambda x:x.getid())
        with open(filename, "w") as savefile:
            #Write the number of nodes in the tree
            savefile.write("%s\n"%(len(nodelist)))

            #Write the basic node data for each node to the file
            for node in nodelist:
                nodeID = node.getid()
                nodescore = node.getscore()
                nodedata = node.getdata()
                datastr = str(nodedata)
                nodeParent = node.getparent()
                if nodeParent is not None:
                    parentID = nodeParent.getid()
                else:
                    parentID = None
                savefile.write("%s %s %s | %s\n"%(nodeID, parentID, nodescore, datastr))

            savefile.write("\n")

            #Write detailed playout information 
            for node in nodelist:
                nodeID = node.getid()
                slist, elist = node.getplayouts()



                usedlist = node.getusedlist()
                nPlayouts = len(slist)
                savefile.write("Node %s %s\n"%(nodeID, nPlayouts))
                for struct, energy in zip(slist, elist):
                    if isinstance(struct,list):
                        outstr = ' '.join([str(x) for x in struct])
                        savefile.write("   %s | %s\n"%(energy, outstr))
                    else:
                        try:
                            savefile.write("   %s | %s\n"%(energy, str(struct)))
                        except:
                            savefile.write("   %s | %s\n"%(energy, struct))

                usedstr = '   ' + ' '.join([str(x) for x in usedlist])
                usedstr = usedstr + "\n"
                savefile.write(usedstr)
                savefile.write("\n")

            #Begin Writing the playout node data for each node to the file

    #-------------------------------------------------------------------------------------------------------------------



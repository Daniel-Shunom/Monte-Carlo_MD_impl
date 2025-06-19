"""
.. module::  Node.py
   :platform: Unix, Windows, Linux
   :synopsis: Node Defintion and Manipulation Opperations

.. moduleauthor:: Troy D Loeffler <tloeffler@anl.gov>

"""

from copy import copy

#================================================
class ObjectiveNotDefined(Exception):
    pass
#================================================
class DataObject(Exception):
    pass

#================================================
class Node(object):
    '''
    Primary Node Object that defines the various Node Manipulation Opperations
    '''
    #--------------------------------------
    def __init__(self, parentNode=None, data=None, objective=None, computescore=True):
        self.__myid = 0
        self.childlist = []
        self.parent = parentNode
        self.libhits = 0
        if self.parent is None:
            self.status = "head"
            self.depth = 0
        else:
            self.status = "child"
            self.depth = self.parent.getdepth() + 1

        self.data = data  #data object such as a string or structural info
        if data is not None and computescore:
            self.objvalue = self.data.computescore() 

        self.minimized = False

        self.playstructs = []
        self.playenergies = []
        self.usedplays = []
    #--------------------------------------
    def __str__(self):
        return "Node %s"%(self.getid())
    #--------------------------------------
    def detailedstring(self):
        '''
         Returns a string containing this node's ID, score, and a data object string.
        '''
        outstr = "Node:%s  Score:%s  Data:%s\n"%(self.getid(), self.objvalue, self.data)
        for energy, struct in zip(self.playenergies, self.playstructs):
            outstr = outstr + "    %s | %s\n"%(energy, struct)
        usedstr = ' '.join([str(x) for x in self.usedplays])
        outstr = outstr + "  Used Playouts: %s"%(usedstr)
        return outstr

    #--------------------------------------
    def spawnchild(self, curid=0):
        '''
         Creates a brand new child node from scratch and appends it to this node's
         child list.
        '''
        print("\n\nSpawining %s from node %s"%(curid+1, self.getid()))
        try:
            newdata = self.data.perturbate(node=self)
        except:
            newdata = self.data.perturbate()
        child = Node(parentNode = self, data = newdata)
#        degenerate = self.degencheck(child)
#        if degenerate:
#            return None
        child.setid(curid+1)
        self.childlist.append(child)
        return child
    #--------------------------------------
    def spawnfromplayout(self, curid=0, choicefunction=None):
        '''
         Creates a brand new child node from this node's playout data
         '''

        # Checks to see if all the playouts have already been used
        if len(self.usedplays) == len(self.playstructs):
            print("Node %s is out of playouts, starting fresh"%(self.getid()))
            child = self.spawnchild(curid=curid)
            return child


        #If no selection function has been given, default to picking the one with the lowest score.
        if choicefunction is None:
            minval = 1e300
            playindex = -1
            for i, score in enumerate(self.playenergies):
                if minval > score:
                    if i not in self.usedplays:
                        playindex = i
                        minval = score
        else:
            playindex = choicefunction(self, self.playenergies, self.playstructures, self.usedplays)
        if playindex < 0:
            return None
        self.usedplays.append(playindex)

        #Redundancy Guard, prevents the same playout from being used to spawn nodes more than once.
        for i, struct in enumerate(self.playstructs):
            if (struct == self.playstructs[playindex]) and (i not in self.usedplays):
                self.usedplays.append(i)
        self.usedplays = sorted(self.usedplays)

        print("\n\nSpawining %s from node %s using playout %s"%(curid+1, self.getid(), playindex))
        newdata = self.data.newdataobject()
        try:
            newdata.setstructure(self.playstructs[playindex])
        except IndexError:
            print(self.playstructs)
            print(self.playenergies)
            print("Struct Length %s"%(len(self.playstructs)))
            print("Energy Length %s"%(len(self.playenergies)))
            print("Given Index %s"%(playindex))
            raise IndexError

        child = Node(parentNode = self, data = newdata)
        child.setscore(self.playenergies[playindex])
#        degenerate = self.degencheck(child)
#        if degenerate:
#            return None
        child.setid(curid+1)
        self.childlist.append(child)
        return child

    #--------------------------------------
    def spawnfromminimize(self, curid=0):
        '''
         Creates a brand new child node from scratch by taking the current node's structure
         and performing a local optimization procedure on it. 
        '''
        print("\n\nSpawining %s from node %s by minizing"%(curid+1, self.getid()))
        newdata, score = self.data.minimize(node=self)
        child = Node(parentNode = self, data = newdata)
        child.setid(curid+1)
        self.childlist.append(child)
        self.minimized = True
        return child, score
    #--------------------------------------
    def spawnfromcrossover(self, partnernode, curid=0):
        '''
         Creates a brand new child node by performing a Genetic Algorithm like
         cross over opperation.
        '''
        print("\n\nSpawining %s from node %s by crossing over with Node %s"%(curid+1, self.getid(), partnernode.getid()))
        newdata = self.data.crossover(node=self, 
                                      pardata=partnernode.getdata(),
                                      parnode=partnernode, 
                                      )
        child = Node(parentNode = self, data = newdata)
        child.setid(curid+1)
        self.childlist.append(child)
        return child

    #--------------------------------------
    def finddeepest(self):
        '''
          Recursively searches for the deepest node in the tree and returns the depth as an integer.
        '''
        depth = self.depth
        if len(self.childlist) < 1:
            return depth
        for child in self.childlist:
            childdepth = child.finddeepest()
            if childdepth > depth:
                depth = childdepth
        return depth
    #--------------------------------------
    def backpropigate(self, nLayer=0, nWins=1, nTries=1):
        '''
        Recursively propagates the current win rate back up the tree. This function is currently depreciated
        in favor of node function calls that gather this node's playouts and all of it's children's playouts.
        '''
        if nLayer == 0:
            if self.parent is None:
                return
            self.parent.backpropigate(nLayer=nLayer+1, nWins=self.wins, nTries=self.sims)
            return

#        print "Node:%s"%(self.getid()), nLayer, nWins, nTries, self.nSumWins, self.nSumSims
        if self.parent is None:
            return
        self.parent.backpropigate(nLayer=nLayer+1, nWins=nWins, nTries=nTries)

    #--------------------------------------
    def evaluate(self):
        '''
          Direct's this node's data object to compute the score of it's structure. Commonly called after the first time
          a node is initialized
        '''
        self.objvalue = self.data.computescore()

    #--------------------------------------
    def simulate(self, playouts=1):
        '''
         Runs a playout simulation 
        '''
        if playouts == 0:
            return
        energies, structures = self.data.runsim(playouts=playouts, node=self)
        self.playstructs = self.playstructs + structures
        self.playenergies = self.playenergies + energies
#        try:
#        except:
#            energies, structures = self.data.runsim(playouts=playouts)
#        self.libhits += libhits

#        self.visits = len(self.playstructs)
    #--------------------------------------
    def associateplayouts(self):
        '''
         Scraps through previously tried data sets to figure out if any of them
         could have potentially been a playout of this node and then associates them
         with this node.
        '''
        energies, structures = self.data.findplayouts(node=self)
        print("Found %s existing points to add as playouts to Node %s"%(len(energies), self.getid()))
        self.playstructs = structures
        self.playenergies = energies

    #--------------------------------------
    def crossover(self, partner):
        '''
         Genetic Algorithm based cross over move between two nodes. New node is spawned as a hybrid of this node's
         data set and it's partner node's data set.
        '''
        newstruct = self.data.crossover(node=self, parnode=partner)
        return newstruct

    #--------------------------------------
    def loaddatastr(self, instr):
        '''
         Used in reloading a tree from a restart file.  Takes the input string and direct this node's data object to convert
         it to a python object and then loads it into the data object.
        '''

        newdata = self.data.convertstr(instr)
        self.data.setstructure(newdata)
    #--------------------------------------
    def strtodata(self, instr):
        '''
        Used to convert an input string from the restart file into a python based data set that
        can be used internally.  Primarily used in restarting a previous search.
        '''
        return self.data.convertstr(instr)
    #--------------------------------------
    def addplayout(self, playstr, playscore):
        '''
         Appends a playout and it's score to this node.
        '''
        if isinstance(playstr, str):
            playstruct = self.strtodata(playstr)
        else:
            playstruct = playstr
        self.playstructs.append(playstruct)
        self.playenergies.append(playscore)
     #--------------------------------------
    def getplayouts(self):
        '''
         Returns this node's playout scores and structures
        '''
        return self.playstructs, self.playenergies    
     #--------------------------------------
    def getallplayouts(self, layer=0):
        '''
         Returns this node's playout scores and structures along with it's children's
        '''
        totalstructs = []
        totalenergies = []
        totalstructs = totalstructs + self.playstructs
        totalenergies = totalenergies + self.playenergies


        if len(self.childlist) < 1:
            return totalstructs, totalenergies

        for child in self.childlist:
            childstructs, childenergies = child.getallplayouts(layer=layer+1)
            totalstructs = totalstructs + childstructs
            totalenergies = totalenergies + childenergies
        return totalstructs, totalenergies
    #--------------------------------------
    def setid(self, myid):
        '''
         Set's this node's ID number
        '''
        self.__myid = myid
     #--------------------------------------
    def getid(self):
        '''
         Get's this node's ID number
        '''
        return self.__myid
     #--------------------------------------
    def isminned(self):
        '''
         Simply returns the "minimized" flag.  This is to prevent a node from being selected
         multiple times for local minimization. 
        '''

        return self.minimized
     #--------------------------------------
    def addchild(self, newnode):
        newnode.setparent(self)
        self.associatechildren(newnode)

    #--------------------------------------
    def getchildren(self, index=None):
        if index is None:
            return self.childlist
        else:
            return self.childlist[index]
    #--------------------------------------
    def removechild(self, index=0):
        '''
         Purges a child node from this node.  Also can be used to delete an entire branch of a tree as the
         specified child node and all of it's children will be deleted.
        '''
        del self.chidlist[index]
    #--------------------------------------
    def associatechildren(self, innode):
        '''
         Adds an already existing node to become a child of this node.
        '''
        self.childlist.append(innode)

    #--------------------------------------
    def dissociatechild(self, removenode):
        self.childlist.remove(removenode)
        print("Node %s has been dissociated from Node %s"%(removenode.getid(), self.getid()))

    #--------------------------------------
    def setparent(self, parent):
        if parent is not None:
            self.status = "child"
            self.parent = parent
        else:
            self.status = "head"
            self.parent = None

    #--------------------------------------
    def getusedlist(self):
        '''
         Gets a list of playouts that have been previously used for node expansion. This is kept to
         prevent the same playing from being continuously used again and again to create new nodes.
        '''
        return self.usedplays
    #--------------------------------------
    def setusedlist(self, usedlist):
        self.usedplays = usedlist
    #--------------------------------------
    def setdata(self, data):
        '''
         Sets this node's data object.
        '''
        self.data = data
    #--------------------------------------
    def getdata(self):
        '''
         Gets this node's data object.
        '''
        return self.data.datatostr()
    #--------------------------------------
    def createdepth(self, depth=0):
        '''
          Automatically assigns node depth from scratch, but can only be called from the headnode.
          This is primarily used to reconstruct node location information after reloading a tree
          from a restart file as this info is not saved. Can also be used if there is a significant
          change in the tree structure. For example if a tree branch is cut and grafted into another
          section of the tree, this can be used to recompute the new depth information.
        '''
        self.depth = depth
        print("Node %s's depth is now: %s"%(self.getid(), self.depth))
        for child in self.childlist:
            child.createdepth(depth=depth+1)

    #--------------------------------------
    def getdepth(self):
        '''
         Get's how deep this node is imbeded in the Tree Structure.  Depth = 0 corresponds to the headnode.
        '''
        return self.depth
    #--------------------------------------
    def changedepth(self, delta = -1, layer=0):
        '''
         Shifts the depth of the node in the event that it's location in the tree is changed.
        '''
        self.depth += delta
        if len(self.childlist) < 1:
            return
        for child in self.childlist:
            sublist = child.changedepth(delta=delta, layer=layer+1)
        return
    #--------------------------------------
    def makeheadnode(self):
        '''
         Changes this node to become the new headnode.  Dissociates itself from all other nodes.
        '''
        self.status = "head"
        self.parent = None
    #--------------------------------------
    def getstatus(self):
        '''
         Returns if this node is a 'head' node or a 'child' node.
        '''
        return self.status
    #--------------------------------------
    def getuniquenessdata(self, nodelist):
        '''
         Gets the finger printing score of this node by calling its data object to compute it. 
        '''
        uniquescore = self.data.getuniqueness(node=self, nodelist=nodelist)
        return uniquescore
    #--------------------------------------
    def getlibuniquenessdata(self):
        try:
            uniquescore = self.data.getlibuniqueness(node=self)
        except:
            uniquescore = self.data.getlibuniqueness()
        return uniquescore
    #--------------------------------------
    def getlibhits(self, layer=0):
        '''
         Returns a value that tells the tree how many times has a library value been used
         by this node. (IE how many times has this node created a data set that is already known) 
         This is used only in special cases.
        '''
        return self.libhits
        if len(self.childlist) < 1:
            return self.libhits
        selfsum = self.libhits
        for child in self.childlist:
            childval = child.getlibhits(layer=layer+1)
            selfsum += childval
        return selfsum
    #--------------------------------------
    def getnodelist(self):
        '''
         Returns a python list containing this node and all of it's child nodes via recursion.
        '''
        if len(self.childlist) < 1:
            return [self]
        totallist = []
        for child in self.childlist:
            sublist = child.getnodelist()
            totallist = totallist + sublist
        return [self] + totallist

    #--------------------------------------
    def getscore(self):
        return self.objvalue
    #--------------------------------------
    def setscore(self, objective):
        self.objvalue = objective

    #--------------------------------------
    def getvisits(self, layer=0):
        '''
         Returns the number of visits to this node and all decendent nodes.
        '''
        visits = len(self.playstructs)
        if len(self.childlist) < 1:
            return visits

        for child in self.childlist:
            childvisits = child.getvisits(layer=layer+1)
            visits += childvisits
        return visits
    #--------------------------------------
    def getparent(self):
        '''
         Returns the parent node of this node. 
        '''
        return self.parent
    #--------------------------------------
    def getlineage(self, layer=0):
        '''
         Returns a list of nodes starting from the node who initiated the call
         tracing the tree back up to the head node. The list will be sorted
         in terms of depth with the deepest node listed first.
        '''
        if self.status == "head":
            return [self]

        parentlineage = self.parent.getlineage(layer=layer+1)
        outlineage = [self]+parentlineage
        if layer == 0:
            outlineage = sorted(outlineage, key=lambda x:x.getdepth(), reverse=True)

        return outlineage

    #--------------------------------------
    def getnewick(self, level=0):
        '''
         Returns a newick tree structure string. Uses recursion to traverse the node tree.
         Used for visualization of the tree using the Ete software.

         For the entire tree this must be called from the head node
        '''

        if len(self.childlist) < 1:
            return str(self.__myid)
        strlist = []
        for child in self.childlist:
            childstr = child.getnewick(level=level+1)
            strlist.append(childstr)
        newstr = ','.join([str(x) for x in strlist])
        newstr = "(" + newstr + ")" + str(self.getid())
        if level == 0:
            newstr = "(" + newstr + ");"
        return newstr

    #--------------------------------------
    def getratiolist(self, layer=0):
        '''
        Previously used to get the wins/tries ratio, but this function is now begining to be
        depreciated. 
        '''
        if len(self.childlist) < 1:
            return [[self.getid(), self.getscore(), len(self.playenergies) ]]
        totallist = []
        for child in self.childlist:
            sublist = child.getratiolist()
            totallist = totallist + sublist
        return [[self.getid(), self.getscore(), len(self.playenergies)]] + totallist
    #--------------------------------------
    def getscorelist(self, layer=0):
        '''
         Returns a score list of this node and all decendent nodes.
        '''
        if len(self.childlist) < 1:
            return [self.getscore()]
        totallist = []
        for child in self.childlist:
            sublist = child.getratiolist()
            totallist = totallist + sublist
        return [self.getscore()] + totallist

    #--------------------------------------
    def getenergylist(self, layer=0):
        '''
         Returns the playout energies of this node and all decendent nodes.
        '''
        if len(self.childlist) < 1:
            return self.playenergies
        totallist = []
        for child in self.childlist:
            sublist = child.getenergylist(layer=layer+1)
            totallist = totallist + sublist
        return self.playenergies + totallist
    #--------------------------------------
    def getbestplayscore(self, sortfunc=min, layer=0):
        '''
         Searches this node's own score and all of it's children nodes for
         the best playout score. Can be given any arbitrary function.
        '''
        try:
            bestscore = sortfunc(self.playenergies)
        except:
            bestscore = None
            print(self.playenergies)
        if len(self.childlist) < 1:
            return bestscore

        for child in self.childlist:
            childbest = child.getbestscore(layer=layer+1)
            if childbest is None:
                continue
            if childbest < bestscore:
                bestscore = childbest
        return bestscore
    #--------------------------------------
    def getbestscore(self, layer=0):
        '''
         Searches this node's own score and all of it's children nodes for
         the best node score.
        '''
        bestscore = self.getscore()
        if len(self.childlist) < 1:
            return bestscore

        for child in self.childlist:
            childbest = child.getbestscore(layer=layer+1)
            if childbest < bestscore:
                bestscore = childbest
        return bestscore



    #--------------------------------------
#================================================

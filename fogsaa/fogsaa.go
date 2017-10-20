package main

import (
	"github.com/emirpasic/gods/trees/binaryheap"

	d "github.com/hivdb/nucamino/data"
	h "github.com/hivdb/nucamino/scorehandler/general"
	a "github.com/hivdb/nucamino/types/amino"
	n "github.com/hivdb/nucamino/types/nucleic"

	"fmt"
	// "github.com/pkg/profile"
)

const posInf = int((^uint(0)) >> 1)
const negInf = -int((^uint(0))>>1) - 1

const (
	ROOT   = 0   // 0b0000
	A0N001 = 0x1 // 0b0001, +  , (+0, +1)
	A0N011 = 0x3 // 0b0011, ++ , (+0, +2)
	A0N111 = 0x7 // 0b0111, +++, (+0, +3)
	A1N000 = 0x8 // 0b1000, ---, (+1, +0)
	A1N001 = 0x9 // 0b1001, --., (+1, +1)
	A1N010 = 0xa // 0b1010, -.-, (+1, +1)
	A1N011 = 0xb // 0b1011, -.., (+1, +2)
	A1N100 = 0xc // 0b1100, .--, (+1, +1)
	A1N101 = 0xd // 0b1101, .-., (+1, +2)
	A1N110 = 0xe // 0b1110, ..-, (+1, +2)
	A1N111 = 0xf // 0b1111, ..., (+1, +3)

	NA1MASK = 0x4
	NA2MASK = 0x2
	NA3MASK = 0x1
)

const nodeTypesLen = 11

var nodeTypes = []int{
	A0N001,
	A0N011,
	A0N111,
	A1N000,
	A1N001,
	A1N010,
	A1N011,
	A1N100,
	A1N101,
	A1N110,
	A1N111,
}

var nodeGaps = []int{
	/*0: */ -1,
	/*0x1: */ 1,
	/*0x2: */ -1,
	/*0x3: */ 2,
	/*0x4: */ -1,
	/*0x5: */ -1,
	/*0x6: */ -1,
	/*0x7: */ 3,
	/*0x8: */ 3,
	/*0x9: */ 2,
	/*0xa: */ 2,
	/*0xb: */ 1,
	/*0xc: */ 2,
	/*0xd: */ 1,
	/*0xe: */ 1,
	/*0xf: */ 0,
}

var nodeNumNAs = []int{
	/*0x0: */ 0,
	/*0x1: */ 1,
	/*0x2: */ 1,
	/*0x3: */ 2,
	/*0x4: */ 1,
	/*0x5: */ 2,
	/*0x6: */ 2,
	/*0x7: */ 3,
}

var nodeGapOpenings = [0x100]int{}
var nodeIndelCodonOpenings = [0x100]bool{}

func init() {
	// opening of insertion gaps
	for _, i := range []int{0, 0x8, 0x9, 0xa, 0xb, 0xc, 0xd, 0xe, 0xf} {
		for _, j := range []int{0x1, 0x3, 0x7} {
			nodeGapOpenings[(i<<4)|j] = 1
		}
		nodeIndelCodonOpenings[(i<<4)|0x7] = true
	}
	// opening of deletion gaps
	for _, i := range []int{0 /* ins */, 0x1, 0x3, 0x7 /* end with NA */, 0x9, 0xb, 0xd, 0xf} {
		for _, j := range []int{0x8, 0x9, 0xb, 0xc, 0xd, 0xe} {
			nodeGapOpenings[(i<<4)|j] = 1
		}
		nodeGapOpenings[(i<<4)|0xa] = 2
		if i != 0x7 {
			// +++--- is not considered as score bonus
			nodeIndelCodonOpenings[(i<<4)|0x8] = true
		}
	}
}

type Node struct {
	posA            int
	posN            int
	nodeType        int
	minFitnessScore int
	maxFitnessScore int
	minSubScore     int
	maxSubScore     int
	minGapScore     int
	maxGapScore     int
	presentScore    int
	isLeaf          bool
	parent          *Node
}

// type sortableNodes [nodeTypesLen]*Node
//
// func (self sortableNodes) Len() int      { return nodeTypesLen }
// func (self sortableNodes) Swap(i, j int) { *self[i], *self[j] = *self[j], *self[i] }
// func (self sortableNodes) Less(i, j int) bool {
// 	if self[i].isExceeded {
// 		return false
// 	} else if self[j].isExceeded {
// 		return true
// 	}
// 	if self[i].maxFitnessScore > self[j].maxFitnessScore {
// 		return true
// 	}
// 	if self[i].maxFitnessScore == self[j].maxFitnessScore &&
// 		self[i].minFitnessScore >= self[j].minFitnessScore {
// 		return true
// 	}
// 	return false
// }

func nodeComparator(a, b interface{}) int {
	n1 := a.(*Node)
	n2 := b.(*Node)
	cmp := n2.maxFitnessScore - n1.maxFitnessScore
	if cmp == 0 {
		cmp = n2.minFitnessScore - n1.minFitnessScore
	}
	return cmp
}

// Calculate present score (prS) for current node
// prS(self) = prS(parent) + score(self)
func (self *Node) calcPresentScore(alignment *Alignment) {
	var (
		aa            a.AminoAcid
		na1, na2, na3 = n.N, n.N, n.N
		posA          = self.posA
		posN          = self.posN
		sh            = alignment.scoreHandler
		nodeType      = self.nodeType
		parentType    = self.parent.nodeType
		presentScore  = self.parent.presentScore
		// +1 -(((0b0111>>3)<<1)-1) means insertion; -1 -(((0b1000>>3)<<1)-1) means deletion
		indelSign         = -((nodeType >> 3) << 1) + 1
		gaps              = nodeGaps[nodeType]
		gapOpenings       = nodeGapOpenings[(parentType<<4)|nodeType]
		indelCodonOpening = nodeIndelCodonOpenings[(parentType<<4)|nodeType]
	)
	// ppp := posA == 0 && (posN == 0 || posN == 3)
	// if ppp {
	// 	fmt.Println(presentScore)
	// }
	if posA == 0 || posN == 0 /* || posA == alignment.aSeqLen || posN == alignment.nSeqLen */ {
		// no penalty/bonus for initial and trailing gaps
		self.presentScore = presentScore
		return
	}

	/* start gap handling */
	presentScore += gaps*sh.GetGapExtensionScore() + gapOpenings*sh.GetGapOpeningScore()
	// fmt.Printf("-x AA%d,NA%d,%04b,gaps=%d,gapOpenings=%d,pS=%d\n", posA, posN, nodeType, gaps, gapOpenings, presentScore)
	if indelCodonOpening {
		presentScore += sh.GetIndelCodonOpeningScore(posA, indelSign)
	}
	if nodeType == A0N111 || nodeType == A1N000 {
		presentScore += sh.GetIndelCodonExtensionScore()
	}
	// if ppp {
	// 	fmt.Printf("%d,%d,%d\n", gaps, gapOpenings, presentScore)
	// }
	/* end gap handling */

	/* start caculating substitution score */
	// if nodeType>>3 == 1 && nodeType != A1N000 {
	if nodeType == A1N111 {
		aa = alignment.getAA(posA)
		if nodeType&NA3MASK == NA3MASK {
			na3 = alignment.getNA(posN)
			posN--
		}
		if nodeType&NA2MASK == NA2MASK {
			na2 = alignment.getNA(posN)
			posN--
		}
		if nodeType&NA1MASK == NA1MASK {
			na1 = alignment.getNA(posN)
		}
		score, cached := sh.GetCachedSubstitutionScore(na1, na2, na3, aa)
		if !cached {
			score = sh.GetSubstitutionScoreNoCache(na1, na2, na3, aa)
		}
		presentScore += score
		// if ppp {
		// 	fmt.Printf("%04b,%s,%s,%d\n", nodeType, a.ToString(aa), n.WriteString([]n.NucleicAcid{na1, na2, na3}), presentScore)
		// }

	}
	// fmt.Printf("__ %04b, %d, gO%d, gs%d\n", nodeType, presentScore, gapOpenings, gaps)
	/* end cacluating substitution score */
	self.presentScore = presentScore
}

// Calculate fitness score (Tmax and Tmin) for current node
// This function assumes that calcPresentScore() has been called
func (self *Node) calcFitnessScore(alignment *Alignment) {
	if self.isLeaf {
		self.maxFitnessScore = self.presentScore
		self.minFitnessScore = self.presentScore
		return
	}
	var (
		windowSize, numGaps int
		indelSign           int
		gapOpened           = false
		sh                  = alignment.scoreHandler
		fitnessMin          = 0
		fitnessMax          = 0
		posAStart           = self.posA + 1
		posNStart           = self.posN + 1
		aRemainNumNAs       = (alignment.aSeqLen + 1 - posAStart) * 3
		nRemainNumNAs       = alignment.nSeqLen + 1 - posNStart
		gO                  = sh.GetGapOpeningScore()
		gE                  = sh.GetGapExtensionScore()
		// iO                  = negInf
		iE    = sh.GetIndelCodonExtensionScore()
		iOMax = negInf
	)
	if aRemainNumNAs > nRemainNumNAs {
		indelSign = -1
		numGaps = aRemainNumNAs - nRemainNumNAs
		// if numGaps > nRemainNumNAs {
		// 	// local alignment ignore trailing gaps
		// 	numGaps = nRemainNumNAs
		// }
		windowSize = nRemainNumNAs / 3
		gapOpened = (self.nodeType>>3 == 1) && self.nodeType%2 == 0
	} else {
		indelSign = 1
		numGaps = nRemainNumNAs - aRemainNumNAs
		// if numGaps > aRemainNumNAs {
		// 	// local alignment ignore trailing gaps
		// 	numGaps = aRemainNumNAs
		// }
		windowSize = aRemainNumNAs / 3
		gapOpened = (self.nodeType>>3 == 0) && self.nodeType%2 == 1
	}

	// add up positional scores
	// if self.nodeType == A1N000 {
	// 	fmt.Printf("%d-%d, %d, %d\n", posAStart, posAEnd, nRemainNumNAs, aRemainNumNAs)
	// }
	self.minSubScore = alignment.getMinSubstitutionScore(posAStart, windowSize, numGaps)
	self.maxSubScore = alignment.getMaxSubstitutionScore(posAStart, windowSize)
	fitnessMin += self.minSubScore
	fitnessMax += self.maxSubScore
	/*for posA := posAStart; posA < posAEnd; posA++ {
		iO = sh.GetIndelCodonOpeningScore(posA, indelSign)
		if iO > iOMax {
			iOMax = iO
		}
	}*/
	// TODO: temporarily
	iOMax = indelSign * 0

	// gap penalties for Tmin
	// This fomular is satisfiable if:
	// 1. gO <= 0
	// 2. gE <= 0
	// 3. iE >= 0
	// 4. gO + iO + 3gE + iE >= 3gO + 3gE
	// TODO: check if these conditions were satisfied
	maxGapsInWindow := windowSize / 2 * 3
	if numGaps < maxGapsInWindow {
		maxGapsInWindow = numGaps
	}
	// fitnessMin += maxGapsInWindow*gO + numGaps*gE
	fitnessMin += maxGapsInWindow * (gO + gE)
	if numGaps > 0 && gapOpened {
		fitnessMin -= gO
	}

	// gap penalties & indel codon bonus for Tmax
	// This fomular is satisfiable if:
	// 1. gO <= 0
	// 2. gE <= 0
	// 3. iE >= 0
	// 4. iOMax <= |gO|
	// TODO: check if these conditions were satisfied
	//fmt.Printf("-1 %d,%d,%d,%d,%d, %04b\n", self.presentScore, fitnessMax, fitnessMin, aRemainNumNAs, nRemainNumNAs, self.nodeType)

	// no score/penalty for trailing gaps
	// fitnessMax += numGaps*gE + numGaps/3*iE
	fitnessMax += 0*iE + 0*iOMax
	// fmt.Printf("== |(%d)-(%d)|=%d\n", nRemainNumNAs, aRemainNumNAs, numGaps)
	// if !gapOpened {
	// 	if numGaps > 0 {
	// 		fitnessMax += gO
	// 	}
	// 	if numGaps > 2 && iOMax > 0 {
	// 		fitnessMax += iOMax
	// 	}
	// }
	// fmt.Printf("%d,%d,%d,%d,%d,%t,%04b\n", self.presentScore, fitnessMax, fitnessMin, aRemainNumNAs, nRemainNumNAs, gapOpened, self.nodeType)
	self.minGapScore = fitnessMin - self.minSubScore
	self.maxGapScore = fitnessMax - self.maxSubScore

	self.maxFitnessScore = self.presentScore + fitnessMax
	self.minFitnessScore = self.presentScore + fitnessMin
}

func (self *Node) expand(alignment *Alignment) [nodeTypesLen]*Node {
	var (
		node     *Node
		children = [nodeTypesLen]*Node{}
	)
	if self.posA == alignment.aSeqLen || self.posN == alignment.nSeqLen {
		// do not expand if reached the ends
		return children
	}
	for i, nodeType := range nodeTypes {
		if self.nodeType != ROOT && self.nodeType>>3 == 0 && nodeType>>2 == 0x2 {
			// **+-** is not allowed
			continue
		} else if self.nodeType>>3 == 0x1 && self.nodeType%2 == 0 && nodeType>>3 == 0 {
			// **-+** is alow not allowed
			continue
		}
		var (
			posA   = self.posA + nodeType>>3
			posN   = self.posN + nodeNumNAs[nodeType&^0x8]
			isLeaf = false
		)
		if posA > alignment.aSeqLen || posN > alignment.nSeqLen {
			// position exceeded
			continue
		} else if posA == alignment.aSeqLen || posN == alignment.nSeqLen {
			isLeaf = true
		}
		children[i] = &Node{
			posA:     posA,
			posN:     posN,
			nodeType: nodeType,
			parent:   self,
			isLeaf:   isLeaf,
		}
		node = children[i]
		node.calcPresentScore(alignment)
		node.calcFitnessScore(alignment)
		// fmt.Printf("- AA%d,NA%d,Type%04b,score%d, [%d,%d], [%d,%d], [%d,%d]\n", node.posA, node.posN, node.nodeType, node.presentScore, node.minFitnessScore, node.maxFitnessScore, node.minSubScore, node.maxSubScore, node.minGapScore, node.maxGapScore)
	}
	// sort.Sort(sortableNodes(children))
	return children
}

type Alignment struct {
	nSeq                    []n.NucleicAcid
	aSeq                    []a.AminoAcid
	nSeqLen                 int
	aSeqLen                 int
	maxSubstitutionScore    []int
	minSubstitutionScore    []int
	maxSubstitutionScoreSum int
	minSubstitutionScoreSum int
	scoreHandler            *h.GeneralScoreHandler
	rootNode                *Node
	optimalLeaf             *Node
}

func New(nSeq []n.NucleicAcid, aSeq []a.AminoAcid, scoreHandler *h.GeneralScoreHandler) *Alignment {
	var (
		nSeqLen                 = len(nSeq)
		aSeqLen                 = len(aSeq)
		maxSubstitutionScore    = make([]int, aSeqLen)
		minSubstitutionScore    = make([]int, aSeqLen)
		maxSubstitutionScoreSum = 0
		minSubstitutionScoreSum = 0
		rootNode                = Node{
			posA:         0,
			posN:         0,
			nodeType:     ROOT,
			presentScore: 0,
		}
	)
	for i, aa := range aSeq {
		score := scoreHandler.GetMaxSubstitutionScore(aa)
		maxSubstitutionScore[i] = score
		maxSubstitutionScoreSum += score
		score = scoreHandler.GetMinSubstitutionScore(aa)
		minSubstitutionScore[i] = score
		minSubstitutionScoreSum += score
	}

	return &Alignment{
		nSeq:                    nSeq,
		aSeq:                    aSeq,
		nSeqLen:                 nSeqLen,
		aSeqLen:                 aSeqLen,
		maxSubstitutionScore:    maxSubstitutionScore,
		minSubstitutionScore:    minSubstitutionScore,
		maxSubstitutionScoreSum: maxSubstitutionScoreSum,
		minSubstitutionScoreSum: minSubstitutionScoreSum,
		scoreHandler:            scoreHandler,
		rootNode:                &rootNode,
	}
}

func (self *Alignment) getMaxSubstitutionScore(startPosA, windowSize int) int {
	var (
		// maxScoreSum int
		leftPtr = startPosA - 1
		// rightPtr    = leftPtr + windowSize
		// scoreSum    = self.maxSubstitutionScoreSum
		// gO          = self.scoreHandler.GetGapOpeningScore()
	)
	// scoreSum += gO
	// for i := 0; i < leftPtr; i++ {
	// 	scoreSum -= self.maxSubstitutionScore[i]
	// }

	if windowSize >= self.aSeqLen-leftPtr {
		windowSize = self.aSeqLen - leftPtr
		// if rightPtr == self.aSeqLen {
		// 	scoreSum -= gO
		// }
		// return scoreSum
	}
	return windowSize * self.scoreHandler.GetMaxSubstitutionScore(a.W)

	/*for i := self.aSeqLen - 1; i >= leftPtr+windowSize; i-- {
		scoreSum -= self.maxSubstitutionScore[i]
	}
	maxScoreSum = scoreSum
	// var l, r int
	for rightPtr < self.aSeqLen {
		scoreSum += self.maxSubstitutionScore[rightPtr] - self.maxSubstitutionScore[leftPtr]
		leftPtr++
		rightPtr++
		if rightPtr == self.aSeqLen {
			scoreSum -= gO
		}
		if scoreSum > maxScoreSum {
			// l = leftPtr
			// r = rightPtr
			maxScoreSum = scoreSum
		}
	}
	// fmt.Printf("= %d %d %d %d %d\n", startPosA, windowSize, maxScoreSum, l, r)
	return maxScoreSum*/

}

func (self *Alignment) getMinSubstitutionScore(startPosA, windowSize, numGaps int) int {
	windowSize -= numGaps
	if windowSize < 0 {
		return 0
	}
	var (
		// minScoreSum int
		leftPtr = startPosA - 1
		// rightPtr    = leftPtr + windowSize
		// scoreSum    = self.minSubstitutionScoreSum
		// gO          = self.scoreHandler.GetGapOpeningScore()
	)
	// scoreSum += gO
	// for i := 0; i < leftPtr; i++ {
	// 	scoreSum -= self.minSubstitutionScore[i]
	// }

	if windowSize >= self.aSeqLen-leftPtr {
		windowSize = self.aSeqLen - leftPtr
		// if rightPtr == self.aSeqLen {
		// 	scoreSum -= gO
		// }
		// return scoreSum
	}
	return windowSize * self.scoreHandler.GetMinSubstitutionScore(a.W)

	/*for i := self.aSeqLen - 1; i >= leftPtr+windowSize; i-- {
		scoreSum -= self.minSubstitutionScore[i]
	}
	minScoreSum = scoreSum
	for rightPtr < self.aSeqLen {
		scoreSum += self.minSubstitutionScore[rightPtr] - self.minSubstitutionScore[leftPtr]
		leftPtr++
		rightPtr++
		// if rightPtr == self.aSeqLen {
		// 	scoreSum -= gO
		// }
		if scoreSum < minScoreSum {
			minScoreSum = scoreSum
		}
	}
	return minScoreSum*/
}

func (self *Alignment) getNA(posN int) n.NucleicAcid {
	return self.nSeq[posN-1]
}

func (self *Alignment) getAA(posA int) a.AminoAcid {
	return self.aSeq[posA-1]
}

var count = 0

/*func (self *Alignment) alignFrom(curNode *Node, depth int) {
	var (
		topChild, nextChild *Node
	)
	// if count > 55 {
	// if curNode.posN != 300 {
	//fmt.Printf("AA%d,NA%d,Type%04b,score%d\n", curNode.posA, curNode.posN, curNode.nodeType, curNode.presentScore)
	// }
	// }
	count++
	if count == 240 {
		panic(0)
	}
	if curNode.isLeaf {
		curNode.optimalScore = curNode.presentScore
		return
	}
	curNode.expand(self)

	// if count > 55 {
	// for i := 0; i < nodeTypesLen; i++ {
	// 	nextChild = curNode.children[i]
	// 	fmt.Printf("- AA%d,NA%d,Type%04b,%d,%d-%d\n", nextChild.posA, nextChild.posN, nextChild.nodeType, nextChild.presentScore, nextChild.maxFitnessScore, nextChild.minFitnessScore)
	// }
	// }

	for i := 0; i < nodeTypesLen; i++ {
		nextChild = curNode.children[i]
		if nextChild.isExceeded {
			continue
		}
		if topChild != nil && topChild.optimalScore > nextChild.maxFitnessScore {
			// prune remaining branches
			fmt.Printf("pruned at %d (%04b), depth %d, %d > %d\n", i, topChild.nodeType, depth, topChild.optimalScore, nextChild.maxFitnessScore)
			break
		}
		self.alignFrom(nextChild, depth+1)
		if topChild == nil || nextChild.optimalScore > topChild.optimalScore {
			topChild = nextChild
		}
	}
	// if depth == 99 {
	// 	fmt.Printf("AA%d,NA%d,Type%04b,score%d\n", topChild.posA, topChild.posN, topChild.nodeType, topChild.optimalScore)
	// }
	curNode.optimalChild = topChild
	curNode.optimalScore = topChild.optimalScore
	fmt.Printf("AA%d,NA%d,Type%04b,score%d\n", curNode.posA, curNode.posN, curNode.nodeType, curNode.presentScore)
}*/

func (self *Alignment) getMatrixIndex(posA int, posN int) int {
	return 2 * (self.aSeqLen*posN + posA)
}

func (self *Alignment) align() {
	var (
		present      = true
		curNode      *Node
		curNodeIf    interface{}
		optimalLeaf  *Node
		queue        = binaryheap.NewWith(nodeComparator)
		bestNodesMap = map[int]*Node{}
	)
	queue.Push(self.rootNode)
	times := 0
	for present {
		curNodeIf, present = queue.Pop()
		if !present {
			break
		}
		curNode = curNodeIf.(*Node)
		mtIdx := self.getMatrixIndex(curNode.posA, curNode.posN)
		bestNode := bestNodesMap[mtIdx]
		if bestNode == nil || nodeComparator(bestNode, curNode) > 0 { // curNode > bestNode
			bestNodesMap[mtIdx] = curNode
		} else {
			continue
		}
		times++
		if curNode.isLeaf {
			if optimalLeaf == nil || optimalLeaf.presentScore < curNode.presentScore {
				optimalLeaf = curNode
			}
		}
		children := curNode.expand(self)
		for _, child := range children {
			if child == nil {
				continue
			}
			mtIdx := self.getMatrixIndex(child.posA, child.posN)
			bestNode := bestNodesMap[mtIdx]
			if bestNode == nil || nodeComparator(bestNode, child) > 0 { // child > bestNode
				// do not replace the original bestNode since the child node hasn't been expanded
				queue.Push(child)
			}
		}
		if optimalLeaf != nil && optimalLeaf.presentScore > curNode.maxFitnessScore {
			break
		}
		// fmt.Printf("AA%d,NA%d,Type%04b,score%d, [%d,%d], [%d,%d], [%d,%d]\n", curNode.posA, curNode.posN, curNode.nodeType, curNode.presentScore, curNode.minFitnessScore, curNode.maxFitnessScore, curNode.minSubScore, curNode.maxSubScore, curNode.minGapScore, curNode.maxGapScore)
	}
	// self.alignFrom(self.rootNode, 0)
	// fmt.Printf("AA%d,NA%d,Type%04b,score%d, [%d,%d], [%d,%d], [%d,%d]\n", optimalLeaf.posA, optimalLeaf.posN, optimalLeaf.nodeType, optimalLeaf.presentScore, optimalLeaf.minFitnessScore, optimalLeaf.maxFitnessScore, optimalLeaf.minSubScore, optimalLeaf.maxSubScore, optimalLeaf.minGapScore, optimalLeaf.maxGapScore)
	print(times, "\n")
	self.optimalLeaf = optimalLeaf
}

func (self *Alignment) GetReport() {
	node := self.optimalLeaf
	for node != self.rootNode {
		fmt.Printf("AA%d,NA%d,Type%04b,Score%d\n", node.posA, node.posN, node.nodeType, node.presentScore)
		node = node.parent
	}

	// var node = self.rootNode
	// for node != nil {
	// 	fmt.Printf("AA%d,NA%d,Type%04b\n", node.posA, node.posN, node.nodeType)
	// }
}

func main() {
	// defer profile.Start(profile.CPUProfile).Stop()
	ref := d.HIV1BRefLookup["POL"]
	seq := n.ReadString("cctcaaatcactctttggcaacgaccccttgtcacagtaaaaataggaggacagctaaaagaagctctattagatacaggagcagatgatacagtattagaagatataaatttgccaggaaaatggaaaccaaaaatgatagggggaattggaggttttatcaaggtaaggcaatatgatcagatacttatagaaatttgtggaaaaaaggctataggtacagtgttagtaggacctacacctgtcaacataattggacgaaatatgttgactcagattggttgtactttaaatttccca")
	// seq = n.ReadString("CCCATTAGTCCTATTGAAACTGTACCAGTAAAATTAAAGCCAGGAATGGATGGCCCAAAAGTTAAACAATGGCCATTGACAGAAGAAAAAATAAAAGCATTAGAAGAGATTTGTGCAGAATTGGAAAAGGAAGGAAAAATTTCAAAAATTGGGCCTGAAAATCCATACAATACTCCAGTATTTGCCATAAAGAAAAAAGACAATAGCAGTGGTAAATGGAGAAAATTAATGGATTTCAGAGAACTTAATAAGAGAACTCAAGATTTCTGGGAAGTTCAATTAGGAATACCACATCCCGCAGGGTTAAAAAAGAAAAAATCAGTAACAGTACTGGATGTGGGTGATGCATATTTTTCAGTTCCCCTAGATGAAGACTTCAGGAAGTATACTGCATTYACCATACCTAGTATAAACAATGAGACACCAGGRATTAGATATCAGTACAATGTGCTYCCACAGGGATGGAAAGGATCACCAGCAATATTCCAAAGTAGCATGACAAARATCTTAGAGCCTTTTAGAAAACAAAATCCAGACATAGTTATCTATCAATACATGGATGATTTGTATGTAGGATCTGACTTAGAAATAGGGCAGCATAGAACAAAAATAGAGGAACTGAGACAACATCTRTGGAAGTGGGGATTTTRCACACCAGACAAAAAACATCAGAAAGAACCTCCATTCCTTTGGATGGGTTATGAACTCCATCCTGATAAATGGACAGTACAGCCTATAGTGCTGCCAGAAAARGACAGYTGGACTGTCAATGACATACAGAAGTTAGTGGGAAAATTGAATTGGGCAAGTCAGATTTATGCAGGGATTAAAGTAAAGCAATTATGTAAACTMCTTAGGGGAACYAAAGCACTAACAGAAGTAGTACCACTAACAGAAGAAGCAGAG")
	alignment := New(seq, ref, h.New(
		/* stopCodonPenalty */ 4,
		/* gapOpeningPenalty */ 10,
		/* gapExtensionPenalty */ 2,
		/* indelCodonOpeningBonus */ 0,
		/* indelCodonExtensionBonus */ 2,
		/* positionalIndelCodonOpeningBonus */ map[int]int{
			// 56prePR + 99PR = 155
			155 + 63: -3,
			155 + 64: -3,
			155 + 65: -9,
			155 + 66: -9,
			155 + 67: -9,
			155 + 68: -3,
			155 + 69: 6, // group all insertions to RT69/POL224
			155 + 70: -3,
			155 + 71: -3,
			155 + 72: -3,
			155 + 73: -3,
		},
		/* isPositionalIndelScoreSupported */ true,
	))
	alignment.align()
	alignment.GetReport()
	fmt.Println("")
}

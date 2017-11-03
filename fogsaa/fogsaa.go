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

	GAP_CLOSED     = 0
	INS_GAP_OPENED = 1
	DEL_GAP_OPENED = 2
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
	// The two fitness scores are only used for ranking nodes
	// but not for pruning any nodes. Instead, the presentScore
	// and the maxPossibleScore are used for pruning.
	posA              int
	posN              int
	gapOpeningStatus  int
	nodeType          int
	lowerFitnessScore int
	upperFitnessScore int
	maxPossibleScore  int
	presentScore      int
	isLeaf            bool
	parent            *Node
}

func getGapOpeningStatus(nodeType int) int {
	// There are 3 gap opening statuses for partial alignment (posA, posN):
	// 1. no gap was opened at the last base; therefore the child must receive gap
	//    opening penalty when starts with "+" (additional NAs) or "-" (missing NAs).
	// 2. a NA insertion gap was opened at the last base; therefore the child must
	//    NOT receive any gap opening penalty when starts with "+".
	// 3. a NA deletion gap was opened at the last base; therefore the child must
	//    NOT receive any gap opening penalty when starts with "-".
	status := GAP_CLOSED
	if nodeType == ROOT {
		return status
	} else if nodeType>>3 == 0 {
		// "insertion"
		status = INS_GAP_OPENED
	} else if nodeType%2 == 0 {
		// "deletion" and end with "-"
		status = DEL_GAP_OPENED
	}
	return status
}

func nodeComparator(a, b interface{}) int {
	n1 := a.(*Node)
	n2 := b.(*Node)
	cmp := n2.upperFitnessScore - n1.upperFitnessScore
	if cmp == 0 {
		cmp = n2.lowerFitnessScore - n1.lowerFitnessScore
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
		self.upperFitnessScore = self.presentScore
		self.lowerFitnessScore = self.presentScore
		self.maxPossibleScore = self.presentScore
		return
	}
	var (
		windowSize, numGaps int
		indelSign           int
		gapOpened           = false
		sh                  = alignment.scoreHandler
		lowerFutureScore    = 0
		upperFutureScore    = 0
		maxFutureScore      = 0
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
	lowerFutureScore += alignment.getMaxNegativeScore(posAStart, windowSize, numGaps)
	upperFutureScore += alignment.getMinPositiveScore(posAStart, windowSize)
	maxFutureScore += alignment.getMaxPositiveScore(posAStart, windowSize)
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
	// lowerFutureScore += maxGapsInWindow*gO + numGaps*gE
	lowerFutureScore += maxGapsInWindow * (gO + gE)
	if numGaps > 0 && gapOpened {
		lowerFutureScore -= gO
	}

	// gap penalties & indel codon bonus for Tmax
	// This fomular is satisfiable if:
	// 1. gO <= 0
	// 2. gE <= 0
	// 3. iE >= 0
	// 4. iOMax <= |gO|
	// TODO: check if these conditions were satisfied
	//fmt.Printf("-1 %d,%d,%d,%d,%d, %04b\n", self.presentScore, upperFutureScore, lowerFutureScore, aRemainNumNAs, nRemainNumNAs, self.nodeType)

	// no score/penalty for trailing gaps
	// upperFutureScore += numGaps*gE + numGaps/3*iE
	upperFutureScore += 0*iE + 0*iOMax
	// fmt.Printf("== |(%d)-(%d)|=%d\n", nRemainNumNAs, aRemainNumNAs, numGaps)
	// if !gapOpened {
	// 	if numGaps > 0 {
	// 		upperFutureScore += gO
	// 	}
	// 	if numGaps > 2 && iOMax > 0 {
	// 		upperFutureScore += iOMax
	// 	}
	// }
	// fmt.Printf("%d,%d,%d,%d,%d,%t,%04b\n", self.presentScore, upperFutureScore, lowerFutureScore, aRemainNumNAs, nRemainNumNAs, gapOpened, self.nodeType)

	self.maxPossibleScore = self.presentScore + maxFutureScore
	self.upperFitnessScore = self.presentScore + upperFutureScore
	self.lowerFitnessScore = self.presentScore + lowerFutureScore
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
		} else if (self.nodeType == A0N001 || self.nodeType == A0N011) && nodeType>>3 == 0 {
			// non-triplet insertion (frameshift) must appear at the end of a gap
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
			posA:             posA,
			posN:             posN,
			gapOpeningStatus: getGapOpeningStatus(nodeType),
			nodeType:         nodeType,
			parent:           self,
			isLeaf:           isLeaf,
		}
		node = children[i]
		node.calcPresentScore(alignment)
		node.calcFitnessScore(alignment)
		// fmt.Printf("- AA%d,NA%d,Type%04b,score%d, [%d,%d], [%d,%d], [%d,%d]\n", node.posA, node.posN, node.nodeType, node.presentScore, node.lowerFitnessScore, node.upperFitnessScore, node.minSubScore, node.maxSubScore, node.minGapScore, node.maxGapScore)
	}
	// sort.Sort(sortableNodes(children))
	return children
}

type Alignment struct {
	nSeq         []n.NucleicAcid
	aSeq         []a.AminoAcid
	nSeqLen      int
	aSeqLen      int
	scoreHandler *h.GeneralScoreHandler
	rootNode     *Node
	optimalLeaf  *Node
}

func New(nSeq []n.NucleicAcid, aSeq []a.AminoAcid, scoreHandler *h.GeneralScoreHandler) *Alignment {
	var (
		nSeqLen  = len(nSeq)
		aSeqLen  = len(aSeq)
		rootNode = Node{
			posA:             0,
			posN:             0,
			gapOpeningStatus: getGapOpeningStatus(ROOT),
			nodeType:         ROOT,
			presentScore:     0,
		}
	)

	return &Alignment{
		nSeq:         nSeq,
		aSeq:         aSeq,
		nSeqLen:      nSeqLen,
		aSeqLen:      aSeqLen,
		scoreHandler: scoreHandler,
		rootNode:     &rootNode,
	}
}

func (self *Alignment) getMinPositiveScore(startPosA, windowSize int) int {
	if windowSize >= self.aSeqLen-startPosA+1 {
		windowSize = self.aSeqLen - startPosA + 1
	}
	return windowSize * self.scoreHandler.GetMinPositiveScore()
}

func (self *Alignment) getMaxPositiveScore(startPosA, windowSize int) int {
	if windowSize >= self.aSeqLen-startPosA+1 {
		windowSize = self.aSeqLen - startPosA + 1
	}
	return windowSize * self.scoreHandler.GetMaxPositiveScore()
}

func (self *Alignment) getMaxNegativeScore(startPosA, windowSize, numGaps int) int {
	windowSize -= numGaps
	if windowSize < 0 {
		return 0
	}
	if windowSize >= self.aSeqLen-startPosA+1 {
		windowSize = self.aSeqLen - startPosA + 1
	}
	return windowSize * self.scoreHandler.GetMaxNegativeScore()
}

func (self *Alignment) getNA(posN int) n.NucleicAcid {
	return self.nSeq[posN-1]
}

func (self *Alignment) getAA(posA int) a.AminoAcid {
	return self.aSeq[posA-1]
}

func (self *Alignment) getNodeIndex(node *Node) int {
	if node.nodeType == ROOT {
		return 0
	}
	return 3*(self.aSeqLen*node.posN+node.posA) + node.gapOpeningStatus
}

func (self *Alignment) align() {
	var (
		present        = true
		curNode        *Node
		curNodeIf      interface{}
		optimalLeaf    *Node
		queue          = binaryheap.NewWith(nodeComparator)
		bestSoFarNodes = map[int]*Node{} // stores expanded best-so-far nodes
	)
	queue.Push(self.rootNode)
	times := 0
	for present {
		curNodeIf, present = queue.Pop()
		if !present {
			break
		}
		curNode = curNodeIf.(*Node)
		nodeIdx := self.getNodeIndex(curNode)

		bestNode := bestSoFarNodes[nodeIdx]
		if bestNode == nil || curNode.presentScore > bestNode.presentScore {
			// curNode is the better path from ROOT to current position & status
			bestSoFarNodes[nodeIdx] = curNode
		} else {
			// PRUNE #1
			// pruned since curNode is no longer promised
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
			nodeIdx := self.getNodeIndex(child)

			bestNode := bestSoFarNodes[nodeIdx]
			// Both present scores are compared to find the better path from ROOT to the
			// child position.
			// The bestNode has same position and gap opening status (see above) as child
			// which indicates their future scores should and will be identical.
			if bestNode == nil || child.presentScore > bestNode.presentScore {
				// do not replace the original bestNode since the child node hasn't been expanded
				queue.Push(child)
			} /* else prune the child since it is worse than expanded node */
			// PRUNE #2
		}
		if optimalLeaf != nil && optimalLeaf.presentScore > curNode.maxPossibleScore {
			// PRUNE #3
			// Unlike "upperFitnessScore", the "maxPossibleScore" takes the "real"
			// maximum positive scores to ensure that no better alignment will be pruned
			// fmt.Printf("AA%d,NA%d,Type%04b,score%d, [%d,%d,%d]\n", curNode.posA, curNode.posN, curNode.nodeType, curNode.presentScore, curNode.lowerFitnessScore, curNode.upperFitnessScore, curNode.maxPossibleScore)
			break
		}
	}
	// fmt.Printf("AA%d,NA%d,Type%04b,score%d, [%d,%d], [%d,%d], [%d,%d]\n", optimalLeaf.posA, optimalLeaf.posN, optimalLeaf.nodeType, optimalLeaf.presentScore, optimalLeaf.lowerFitnessScore, optimalLeaf.upperFitnessScore, optimalLeaf.minSubScore, optimalLeaf.maxSubScore, optimalLeaf.minGapScore, optimalLeaf.maxGapScore)
	print(times, "\n")
	self.optimalLeaf = optimalLeaf
}

func (self *Alignment) GetReport() {
	node := self.optimalLeaf
	for node != self.rootNode {
		fmt.Printf("AA%d,NA%d,Type%04b,Score%d\n", node.posA, node.posN, node.nodeType, node.presentScore)
		node = node.parent
	}
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

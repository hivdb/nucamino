package main

import (
	"bufio"
	"fmt"
	"log"
	"os"
	"strings"
	"sync"
)

import (
	"github.com/hivdb/nucamino/alignment"
	//"github.com/hivdb/nucamino/scorehandler/blosum62"
	d "github.com/hivdb/nucamino/data"
	s "github.com/hivdb/nucamino/scorehandler"
	"github.com/hivdb/nucamino/scorehandler/hiv1b"
	a "github.com/hivdb/nucamino/types/amino"
	n "github.com/hivdb/nucamino/types/nucleic"
	//"github.com/pkg/profile"
)

type Sequence struct {
	Name     string
	Sequence []n.NucleicAcid
}

//const GBFILE = "GB.local.seqs.txt"
const GBFILE = "GB.local.seqs.hiv1.txt"

//const RESULTFILES = "Result.GB.local.seqs.-.%s.tsv"
const RESULTFILES = "Result.GB.local.seqs.hiv1.-.%s.tsv"

const COUNT = 120000

func iterSequences(fileName string, first int) chan Sequence {
	resultChan := make(chan Sequence)
	go func() {
		file, err := os.Open(fileName)
		if err != nil {
			log.Fatal(err)
			return
		}
		scanner := bufio.NewScanner(file)
		for scanner.Scan() {
			line := scanner.Text()
			splitted := strings.SplitN(line, "||", 2)
			name, seqText := splitted[0], splitted[1]
			resultChan <- Sequence{name, n.ReadString(seqText)}
			first--
			if first == 0 {
				break
			}
		}
		close(resultChan)
	}()
	return resultChan
}

/*func test(scoreHandler s.ScoreHandler) *alignment.Alignment {
	cmp := alignment.NewAlignment(
		scoreHandler,
	)
	cmp.CalcScore()
	return cmp
}*/

func getReport(seqName string, r alignment.AlignmentReport) string {
	v := fmt.Sprintf("%s\t%d\t%d\t%d\t%d\t%s\t%s\n",
		seqName,
		r.FirstAA, r.LastAA,
		r.FirstNA, r.LastNA,
		func() string {
			muts := ""
			for _, mut := range r.Mutations {
				muts += mut.ToString() + ","
			}
			return muts
		}(),
		func() string {
			fss := ""
			for _, fs := range r.FrameShifts {
				fss += fs.ToString() + ","
			}
			return fss
		}())
	return v
}

func writeReport(seqName string, r alignment.AlignmentReport, fp *os.File) {
	fp.WriteString(getReport(seqName, r))
}

func main() {
	// defer profile.Start(profile.CPUProfile).Stop()
	var (
		wg         = sync.WaitGroup{}
		threads    = 4
		count      = COUNT
		seqChan    = iterSequences(GBFILE, count)
		seqChan2   = iterSequences(GBFILE, count)
		resultChan = make(chan [3]string)
		resultMap  = make(map[string][3]string)
		genes      = [3]hiv1b.Gene{hiv1b.PR, hiv1b.RT, hiv1b.IN}
		refs       = [3][]a.AminoAcid{d.HIV1BSEQ_PR, d.HIV1BSEQ_RT, d.HIV1BSEQ_IN}
		fps        = [3]*os.File{}
	)
	fps[0], _ = os.Create(fmt.Sprintf(RESULTFILES, "PR"))
	fps[1], _ = os.Create(fmt.Sprintf(RESULTFILES, "RT"))
	fps[2], _ = os.Create(fmt.Sprintf(RESULTFILES, "IN"))
	for i := 0; i < threads; i++ {
		wg.Add(1)
		go func(idx int, rChan chan<- [3]string) {
			var (
				scoreHandlers [3]s.ScoreHandler
			)
			for i, gene := range genes {
				scoreHandlers[i] = hiv1b.NewAsScoreHandler(
					/* gene                		*/ gene,
					/* indelCodonOpeningBonus   */ 4,
					/* indelCodonExtensionBonus */ 0,
					/* stopCodonPenalty    		*/ 4,
					/* gapOpenPenalty      		*/ 10,
					/* gapExtensionPenalty 		*/ 2,
				)
			}
			for seq, ok := <-seqChan; ok; seq, ok = <-seqChan {
				result := [3]string{}
				for i := 0; i < 3; i++ {
					aligned, success := alignment.NewAlignment(seq.Name, seq.Sequence, refs[i], scoreHandlers[i])
					if !success {
						continue
					}
					r := aligned.GetReport()
					result[i] = getReport(seq.Name, r)
				}
				rChan <- result
			}
			wg.Done()
		}(i, resultChan)
	}
	go func(rChan chan<- [3]string) {
		wg.Wait()
		close(rChan)
	}(resultChan)
	for result := range resultChan {
		name := strings.SplitN(result[0], "\t", 2)
		resultMap[name[0]] = result
	}
	for seq := range seqChan2 {
		result := resultMap[seq.Name]
		for i, _ := range genes {
			fps[i].WriteString(result[i])
		}
	}
}

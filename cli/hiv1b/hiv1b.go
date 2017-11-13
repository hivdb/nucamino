package cli

import (
	"bytes"
	json "encoding/json"
	"fmt"
	"github.com/hivdb/nucamino/alignment"
	d "github.com/hivdb/nucamino/data"
	h "github.com/hivdb/nucamino/scorehandler/general"
	a "github.com/hivdb/nucamino/types/amino"
	"github.com/hivdb/nucamino/utils/fastareader"
	"log"
	"os"
	"runtime"
	"sync"
)

type Gene uint8

const (
	GAG Gene = iota
	POL
	GP41
)

var GeneLookup = map[string]Gene{
	"GAG":  GAG,
	"POL":  POL,
	"GP41": GP41,
}

var AllPositionalIndelScores = map[Gene]map[int]int{
	POL: map[int]int{
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
}

type tAlignmentResult struct {
	Name   string
	Report *alignment.AlignmentReport
}

func writeTSV(
	file *os.File, textGenes []string,
	seqs []fastareader.Sequence, resultMap map[string][]tAlignmentResult) {

	genesCount := len(textGenes)
	file.WriteString("Sequence Name")
	for _, textGene := range textGenes {
		file.WriteString("\t" + textGene + " FirstAA")
		file.WriteString("\t" + textGene + " LastAA")
		file.WriteString("\t" + textGene + " FirstNA")
		file.WriteString("\t" + textGene + " LastNA")
		file.WriteString("\t" + textGene + " Mutations")
		file.WriteString("\t" + textGene + " FrameShifts")
	}
	file.WriteString("\n")

	for _, seq := range seqs {
		result := resultMap[seq.Name]
		if result == nil {
			continue
		}
		file.WriteString(seq.Name)
		for i := 0; i < genesCount; i++ {
			r := result[i].Report
			file.WriteString(fmt.Sprintf(
				"\t%d\t%d\t%d\t%d\t%s\t%s",
				r.FirstAA, r.LastAA,
				r.FirstNA, r.LastNA,
				func() string {
					var muts bytes.Buffer
					for _, mut := range r.Mutations {
						muts.WriteString(mut.ToString())
						muts.WriteString(",")
					}
					if muts.Len() > 0 {
						muts.Truncate(muts.Len() - 1)
					}
					return muts.String()
				}(),
				func() string {
					var fss bytes.Buffer
					for _, fs := range r.FrameShifts {
						fss.WriteString(fs.ToString())
						fss.WriteString(",")
					}
					if fss.Len() > 0 {
						fss.Truncate(fss.Len() - 1)
					}
					return fss.String()
				}(),
			))
		}
		file.WriteString("\n")
	}
}

func writeJSON(
	file *os.File, textGenes []string,
	seqs []fastareader.Sequence, resultMap map[string][]tAlignmentResult) {

	finalResultMap := make(map[string][]tAlignmentResult)
	genesCount := len(textGenes)

	for i := 0; i < genesCount; i++ {
		textGene := textGenes[i]
		for _, seq := range seqs {
			seqResult := resultMap[seq.Name]
			if seqResult != nil {
				seqGeneResult := seqResult[i]
				finalResultMap[textGene] = append(finalResultMap[textGene], seqGeneResult)
			}
		}
	}
	result, err := json.MarshalIndent(finalResultMap, "", "  ")
	if err != nil {
		log.Fatal(err)
	}
	file.Write(result)
}

func seqSlice2Chan(s []fastareader.Sequence, bufferSize int) chan fastareader.Sequence {
	c := make(chan fastareader.Sequence, bufferSize)
	go func() {
		for _, item := range s {
			c <- item
		}
		close(c)
	}()
	return c
}

func PerformAlignment(
	inputFileName string,
	outputFileName string,
	textGenes []string,
	indelCodonOpeningBonus int,
	indelCodonExtensionBonus int,
	stopCodonPenalty int,
	gapOpeningPenalty int,
	gapExtensionPenalty int,
	goroutines int,
	outputFormat string, quiet bool) {

	runtime.LockOSThread()
	numCPU := runtime.NumCPU()
	logger := log.New(os.Stderr, "", 0)
	if goroutines == 0 {
		goroutines = numCPU
	}
	if !quiet {
		logger.Printf(
			"%d CPUs were detected. %d goroutines will be created.\n",
			numCPU, goroutines)
	}

	var input, output *os.File
	var err error
	if inputFileName == "-" {
		input = os.Stdin
	} else {
		input, err = os.Open(inputFileName)
		if err != nil {
			logger.Fatal(err)
			return
		}
	}
	if outputFileName == "-" {
		output = os.Stdout
	} else {
		output, err = os.Create(outputFileName)
		if err != nil {
			logger.Fatal(err)
			return
		}
	}

	genesCount := len(textGenes)
	genes := make([]Gene, genesCount)
	refs := make([][]a.AminoAcid, genesCount)
	for i, textGene := range textGenes {
		genes[i] = GeneLookup[textGene]
		refs[i] = d.HIV1BRefLookup[textGene]
	}

	var (
		wg         = sync.WaitGroup{}
		seqs       = fastareader.ReadSequences(input)
		resultChan = make(chan []tAlignmentResult)
		resultMap  = make(map[string][]tAlignmentResult)
	)
	if !quiet {
		logger.Printf("%d sequences were found from the input file.\n", len(seqs))
	}
	var seqChan = seqSlice2Chan(seqs, goroutines*4)
	for i := 0; i < goroutines; i++ {
		wg.Add(1)
		go func(idx int, rChan chan<- []tAlignmentResult) {
			scoreHandlers := make([]*h.GeneralScoreHandler, genesCount)
			for i, gene := range genes {
				positionalIndelScores, isPositionalIndelScoreSupported := AllPositionalIndelScores[gene]
				scoreHandlers[i] = h.New(
					stopCodonPenalty,
					gapOpeningPenalty,
					gapExtensionPenalty,
					indelCodonOpeningBonus,
					indelCodonExtensionBonus,
					positionalIndelScores,
					isPositionalIndelScoreSupported)
			}
			for seq := range seqChan {
				isSimpleAlignment := true
				result := make([]tAlignmentResult, genesCount)
				for i := 0; i < genesCount; i++ {
					aligned := alignment.NewAlignment(seq.Sequence, refs[i], scoreHandlers[i])
					r := aligned.GetReport()
					result[i] = tAlignmentResult{seq.Name, r}
					isSimpleAlignment = isSimpleAlignment && r.IsSimpleAlignment
				}
				rChan <- result
				if !quiet {
					if isSimpleAlignment {
						fmt.Fprintf(os.Stderr, ":")
					} else {
						fmt.Fprintf(os.Stderr, ".")
					}
				}
			}
			wg.Done()
		}(i, resultChan)
	}
	go func(rChan chan<- []tAlignmentResult) {
		wg.Wait()
		if !quiet {
			logger.Printf("\n")
		}
		close(rChan)
	}(resultChan)
	for result := range resultChan {
		resultMap[result[0].Name] = result
	}
	switch outputFormat {
	case "tsv":
		writeTSV(output, textGenes, seqs, resultMap)
		break
	case "json":
		writeJSON(output, textGenes, seqs, resultMap)
		break
	}
	if !quiet && outputFileName != "-" {
		logger.Printf("Created alignment result file %s.", outputFileName)
	}
}

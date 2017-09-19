package cli

import (
	"bytes"
	"fmt"
	"github.com/hivdb/nucamino/alignment"
	d "github.com/hivdb/nucamino/data"
	s "github.com/hivdb/nucamino/scorehandler"
	hiv1bhandler "github.com/hivdb/nucamino/scorehandler/hiv1b"
	a "github.com/hivdb/nucamino/types/amino"
	"github.com/hivdb/nucamino/utils/fastareader"
	"log"
	"os"
	"runtime"
	"sync"
)

func reportToTSVRow(seqName string, r alignment.AlignmentReport) string {
	v := fmt.Sprintf(
		"%s\t%d\t%d\t%d\t%d\t%s\t%s",
		seqName,
		r.FirstAA, r.LastAA,
		r.FirstNA, r.LastNA,
		func() string {
			var muts bytes.Buffer
			for _, mut := range r.Mutations {
				muts.WriteString(mut.ToString())
				muts.WriteString(",")
			}
			muts.Truncate(muts.Len() - 1)
			return muts.String()
		}(),
		func() string {
			var fss bytes.Buffer
			for _, fs := range r.FrameShifts {
				fss.WriteString(fs.ToString())
				fss.WriteString(",")
			}
			fss.Truncate(fss.Len() - 1)
			return fss.String()
		}())
	return v
}

type alignmentResult struct {
	Name   string
	Report *alignment.AlignmentReport
}

func writeTSV(
	file *os.File, textGenes []string,
	seqs []fastareader.Sequence, resultMap map[string][]alignmentResult) {

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
	goroutines int, quiet bool) {

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
	genes := make([]hiv1bhandler.Gene, genesCount)
	refs := make([][]a.AminoAcid, genesCount)
	for i, textGene := range textGenes {
		genes[i] = hiv1bhandler.GeneLookup[textGene]
		refs[i] = d.HIV1BRefLookup[textGene]
	}

	var (
		wg         = sync.WaitGroup{}
		seqs       = fastareader.ReadSequences(input)
		resultChan = make(chan []alignmentResult)
		resultMap  = make(map[string][]alignmentResult)
	)
	if !quiet {
		logger.Printf("%d sequences were found from the input file.\n", len(seqs))
	}
	var seqChan = seqSlice2Chan(seqs, goroutines*4)
	for i := 0; i < goroutines; i++ {
		wg.Add(1)
		go func(idx int, rChan chan<- []alignmentResult) {
			scoreHandlers := make([]s.ScoreHandler, genesCount)
			for i, gene := range genes {
				scoreHandlers[i] = hiv1bhandler.NewAsScoreHandler(
					gene,
					indelCodonOpeningBonus,
					indelCodonExtensionBonus,
					stopCodonPenalty,
					gapOpeningPenalty,
					gapExtensionPenalty)
			}
			for seq := range seqChan {
				result := make([]alignmentResult, genesCount)
				for i := 0; i < genesCount; i++ {
					aligned, success := alignment.NewAlignment(seq.Sequence, refs[i], scoreHandlers[i])
					if !success {
						continue
					}
					r := aligned.GetReport()
					result[i] = alignmentResult{seq.Name, &r}
				}
				rChan <- result
				if !quiet {
					fmt.Fprintf(os.Stderr, ".")
				}
			}
			wg.Done()
		}(i, resultChan)
	}
	go func(rChan chan<- []alignmentResult) {
		wg.Wait()
		if !quiet {
			logger.Printf("\n")
		}
		close(rChan)
	}(resultChan)
	for result := range resultChan {
		resultMap[result[0].Name] = result
	}
	writeTSV(output, textGenes, seqs, resultMap)
	if !quiet && outputFileName != "-" {
		logger.Printf("Created alignment result file %s.", outputFileName)
	}
}

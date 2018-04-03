package cli

import (
	"bytes"
	"encoding/json"
	"fmt"
	"github.com/hivdb/nucamino/alignment"
	ap "github.com/hivdb/nucamino/alignmentprofile"
	h "github.com/hivdb/nucamino/scorehandler/general"
	a "github.com/hivdb/nucamino/types/amino"
	"github.com/hivdb/nucamino/utils/fastareader"
	"log"
	"os"
	"runtime"
	"sync"
)

type AlignmentResult struct {
	Name   string
	Report *alignment.AlignmentReport
	Error  string
	Err    error
}

func validOutputFormat(format string) bool {
	validFormats := []string{"json", "tsv"}
	for _, validFormat := range validFormats {
		if format == validFormat {
			return true
		}
	}
	return false
}

func writeTSV(
	file *os.File, textGenes []string,
	seqs []fastareader.Sequence, resultMap map[string][]AlignmentResult) {

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
			err := result[i].Err
			if err != nil {
				file.WriteString("\tNA\tNA\tNA\tNA\tNA\tNA")
				continue
			}
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
	seqs []fastareader.Sequence, resultMap map[string][]AlignmentResult) {

	finalResultMap := make(map[string][]AlignmentResult)
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
	outputFormat string,
	textGenes []string,
	goroutines int,
	quiet bool,
	alignmentProfile ap.AlignmentProfile) error {

	// Check output format
	if !validOutputFormat(outputFormat) {
		err := fmt.Errorf("Unknown output format %v. Options are: tsv, json", outputFormat)
		return err
	}

	// Configure runtime
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

	// Prepare input and output files
	var input, output *os.File
	var err error

	if inputFileName == "-" {
		input = os.Stdin
	} else {
		input, err = os.Open(inputFileName)
		if err != nil {
			return err
		}
	}
	if outputFileName == "-" {
		output = os.Stdout
	} else {
		output, err = os.Create(outputFileName)
		if err != nil {
			return err
		}
	}

	genesCount := len(textGenes)
	genes := make([]ap.Gene, genesCount)
	refs := make([][]a.AminoAcid, genesCount)
	for i, textGene := range textGenes {
		genes[i] = ap.Gene(textGene)
		refs[i] = alignmentProfile.ReferenceSequences[genes[i]]
	}

	var (
		wg         = sync.WaitGroup{}
		seqs       = fastareader.ReadSequences(input)
		resultChan = make(chan []AlignmentResult)
		resultMap  = make(map[string][]AlignmentResult)
	)
	if !quiet {
		logger.Printf("%d sequences were found from the input file.\n", len(seqs))
	}

	var seqChan = seqSlice2Chan(seqs, goroutines*4)
	for i := 0; i < goroutines; i++ {
		wg.Add(1)
		go func(idx int, rChan chan<- []AlignmentResult) {
			scoreHandlers := make([]*h.GeneralScoreHandler, genesCount)
			for i, gene := range genes {
				scoreHandlers[i] = h.New(gene, alignmentProfile)
			}
			for seq := range seqChan {
				isSimpleAlignment := true
				result := make([]AlignmentResult, genesCount)
				for i := 0; i < genesCount; i++ {
					aligned, err := alignment.NewAlignment(seq.Sequence, refs[i], scoreHandlers[i])
					if err != nil {
						result[i] = AlignmentResult{seq.Name, nil, err.Error(), err}
					} else {
						r := aligned.GetReport()
						result[i] = AlignmentResult{seq.Name, r, "", nil}
						isSimpleAlignment = isSimpleAlignment && r.IsSimpleAlignment
					}
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
	go func(rChan chan<- []AlignmentResult) {
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
	return nil
}

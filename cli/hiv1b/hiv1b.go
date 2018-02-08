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

var AllPositionalIndelScores = map[Gene]map[int][2]int{
	GAG: map[int][2]int{
		111: [2]int{-5, 0},
		112: [2]int{-5, 0},
		113: [2]int{11, 0},
		114: [2]int{-5, 0},

		115: [2]int{-6, 0},
		116: [2]int{-6, 0},
		117: [2]int{15, 0},
		118: [2]int{-6, 0},
		119: [2]int{-6, 0},

		124: [2]int{-6, 0},
		125: [2]int{-6, 0},
		126: [2]int{15, 0},
		127: [2]int{-6, 0},

		// MA/CA
		128: [2]int{-6, -2},
		129: [2]int{-2, -2},
		130: [2]int{-2, -2},
		131: [2]int{-2, -2},
		132: [2]int{-2, -2},
		133: [2]int{-2, -2},
		134: [2]int{-2, -2},
		135: [2]int{-2, -2},
		136: [2]int{-2, -2},
		137: [2]int{-2, -2},

		249: [2]int{-6, 0},
		250: [2]int{-6, 0},
		251: [2]int{15, 0},
		252: [2]int{-6, 0},
		253: [2]int{-6, 0},

		// CA/SP1
		359: [2]int{-2, -2},
		360: [2]int{-2, -2},
		361: [2]int{-2, -2},
		362: [2]int{-2, -2},
		363: [2]int{-2, -2},
		364: [2]int{-2, -2},
		365: [2]int{-2, -2},
		366: [2]int{-5, -2},
		367: [2]int{-5, -2},
		368: [2]int{11, -2},

		369: [2]int{-6, 0},
		370: [2]int{-6, 0},
		371: [2]int{15, 0},
		372: [2]int{-6, 0},

		// SP1/NC
		373: [2]int{-6, -2},
		374: [2]int{-2, -2},
		375: [2]int{-2, -2},
		376: [2]int{-2, -2},
		377: [2]int{-2, -2},
		378: [2]int{-2, -2},
		379: [2]int{-2, -2},
		380: [2]int{-2, -2},
		381: [2]int{-5, -2},
		382: [2]int{-5, -2},

		383: [2]int{11, 0},
		384: [2]int{-5, 0},
		385: [2]int{-5, 0},
		386: [2]int{-3, 0},

		390: [2]int{0, 1},

		424: [2]int{-6, 0},
		425: [2]int{-6, 0},
		426: [2]int{14, 0},
		427: [2]int{-6, 0},

		// NC/SP2
		428: [2]int{-6, -2},
		429: [2]int{-2, -2},
		430: [2]int{-2, -2},
		431: [2]int{-2, -2},
		432: [2]int{-2, -2},
		433: [2]int{-2, -2},
		434: [2]int{-2, -2},
		435: [2]int{-2, -2},
		436: [2]int{-2, -2},
		437: [2]int{-2, -2},

		438: [2]int{-2, 0},
		439: [2]int{-4, 0},
		440: [2]int{-4, 0},
		441: [2]int{9, 0},
		442: [2]int{-4, 0},
		443: [2]int{-4, 0},

		// SP2/p6
		444: [2]int{-2, -2},
		445: [2]int{-2, -2},
		446: [2]int{-2, -2},
		447: [2]int{-2, -2},
		448: [2]int{-2, -2},
		449: [2]int{-2, -2},
		450: [2]int{-2, -2},
		451: [2]int{-5, -2},
		452: [2]int{-5, -2},
		453: [2]int{11, -2},

		454: [2]int{-5, 0},
		455: [2]int{-6, 0},
		456: [2]int{14, 0},
		457: [2]int{-6, 0},
		458: [2]int{-5, 0},

		465: [2]int{0, 1},
		466: [2]int{0, 1},

		467: [2]int{-3, 0},
		468: [2]int{-3, 0},
		469: [2]int{-3, 0},
		470: [2]int{-3, 0},
		471: [2]int{-4, 0},
		472: [2]int{-5, 0},
		473: [2]int{14, 0},
		474: [2]int{-6, 0},
		475: [2]int{-6, 0},
		476: [2]int{-5, 0},
		477: [2]int{-4, 0},

		478: [2]int{-4, 0},
		479: [2]int{-5, 0},
		480: [2]int{-5, 0},
		481: [2]int{-6, 0},
		482: [2]int{14, 0},

		// p6/PR
		483: [2]int{-5, -2},
		484: [2]int{-5, -2},
		485: [2]int{-2, -2},
		486: [2]int{-2, -2},
		487: [2]int{-2, -2},
		488: [2]int{-2, -2},
		489: [2]int{-2, -2},
		490: [2]int{-2, -2},
		491: [2]int{-2, -2},
		492: [2]int{-2, -2},
		493: [2]int{-2, -2},
	},
	POL: map[int][2]int{
		// 56prePR + 99PR = 155
		155 + 63:  [2]int{-5, 0},
		-155 - 63: [2]int{-5, 0},
		155 + 64:  [2]int{-5, 0},
		-155 - 64: [2]int{-5, 0},
		155 + 65:  [2]int{-7, 0},
		155 + 66:  [2]int{-7, 0},
		155 + 67:  [2]int{-7, 0},
		155 + 68:  [2]int{-3, 0},
		155 + 69:  [2]int{18, -3}, // group all insertions to RT69/POL224
		155 + 70:  [2]int{-3, 0},
		155 + 71:  [2]int{-3, 0},
		155 + 72:  [2]int{-3, 0},
		155 + 73:  [2]int{-3, 0},
	},
}

type tAlignmentResult struct {
	Name   string
	Report *alignment.AlignmentReport
	Error  string
	err    error
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
			err := result[i].err
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
					aligned, err := alignment.NewAlignment(seq.Sequence, refs[i], scoreHandlers[i])
					if err != nil {
						result[i] = tAlignmentResult{seq.Name, nil, err.Error(), err}
					} else {
						r := aligned.GetReport()
						result[i] = tAlignmentResult{seq.Name, r, "", nil}
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

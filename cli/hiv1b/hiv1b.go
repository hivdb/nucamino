package hiv1b

import (
	"fmt"
	"github.com/hivdb/nucamino/alignment"
	cli "github.com/hivdb/nucamino/cli/cli"
	d "github.com/hivdb/nucamino/data"
	h "github.com/hivdb/nucamino/scorehandler/general"
	a "github.com/hivdb/nucamino/types/amino"
	"github.com/hivdb/nucamino/utils/fastareader"
	"log"
	"os"
	"runtime"
	"sync"
)

var AllPositionalIndelScores = map[cli.Gene]map[int][2]int{
	cli.GAG: map[int][2]int{
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
	cli.POL: map[int][2]int{
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
	genes := make([]cli.Gene, genesCount)
	refs := make([][]a.AminoAcid, genesCount)
	for i, textGene := range textGenes {
		genes[i] = cli.GeneLookup[textGene]
		refs[i] = d.HIV1BRefLookup[textGene]
	}

	var (
		wg         = sync.WaitGroup{}
		seqs       = fastareader.ReadSequences(input)
		resultChan = make(chan []cli.AlignmentResult)
		resultMap  = make(map[string][]cli.AlignmentResult)
	)
	if !quiet {
		logger.Printf("%d sequences were found from the input file.\n", len(seqs))
	}
	var seqChan = cli.SeqSlice2Chan(seqs, goroutines*4)
	for i := 0; i < goroutines; i++ {
		wg.Add(1)
		go func(idx int, rChan chan<- []cli.AlignmentResult) {
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
				result := make([]cli.AlignmentResult, genesCount)
				for i := 0; i < genesCount; i++ {
					aligned, err := alignment.NewAlignment(seq.Sequence, refs[i], scoreHandlers[i])
					if err != nil {
						result[i] = cli.AlignmentResult{seq.Name, nil, err.Error(), err}
					} else {
						r := aligned.GetReport()
						result[i] = cli.AlignmentResult{seq.Name, r, "", nil}
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
	go func(rChan chan<- []cli.AlignmentResult) {
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
		cli.WriteTSV(output, textGenes, seqs, resultMap)
		break
	case "json":
		cli.WriteJSON(output, textGenes, seqs, resultMap)
		break
	}
	if !quiet && outputFileName != "-" {
		logger.Printf("Created alignment result file %s.", outputFileName)
	}
}

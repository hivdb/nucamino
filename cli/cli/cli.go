package cli

import (
	"bytes"
	"encoding/json"
	"fmt"
	"github.com/hivdb/nucamino/alignment"
	"github.com/hivdb/nucamino/utils/fastareader"
	"log"
	"os"
)

type Gene uint8

const (
	// hiv1b genes
	GAG Gene = iota
	POL
	GP41
)

var GeneLookup = map[string]Gene{
	"GAG":  GAG,
	"POL":  POL,
	"GP41": GP41,
}


type AlignmentResult struct {
	Name   string
	Report *alignment.AlignmentReport
	Error  string
	Err    error
}


func WriteTSV(
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

func WriteJSON(
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

func SeqSlice2Chan(s []fastareader.Sequence, bufferSize int) chan fastareader.Sequence {
	c := make(chan fastareader.Sequence, bufferSize)
	go func() {
		for _, item := range s {
			c <- item
		}
		close(c)
	}()
	return c
}

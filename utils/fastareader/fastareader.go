package fastareader

import (
	"bufio"
	"bytes"
	"fmt"
	n "github.com/hivdb/nucamino/types/nucleic"
	"io"
	"strings"
)

type Sequence struct {
	Name     string
	Sequence []n.NucleicAcid
}

func makeSequence(name string, seqText string) Sequence {
	return Sequence{name, n.ReadString(seqText)}
}

func ReadSequences(reader io.Reader) []Sequence {
	results := make([]Sequence, 0, 20)
	name := ""
	var seqBuffer bytes.Buffer
	seqCount := 0
	scanner := bufio.NewScanner(reader)
	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, ";") || strings.HasPrefix(line, "#") {
			continue
		} else if strings.HasPrefix(line, ">") {
			if name != "" {
				results = append(
					results, makeSequence(name, seqBuffer.String()))
				seqBuffer.Reset()
			}
			seqCount++
			name = strings.TrimSpace(strings.TrimPrefix(line, ">"))
			if name == "" {
				name = fmt.Sprintf("unnamed sequence %d", seqCount)
			}
		} else {
			seqBuffer.WriteString(strings.TrimSpace(line))
		}
	}
	if name != "" || seqBuffer.Len() > 0 {
		if name == "" {
			name = "unnamed sequence"
		}
		results = append(results, makeSequence(name, seqBuffer.String()))
	}
	return results
}

package main

import "os"
import "github.com/jessevdk/go-flags"

type FileOptions struct {
	Input  flags.Filename `short:"i" long:"input" required:"yes" value-name:"INPUT" description:"FASTA file contains one or more DNA sequences" default:"-"`
	Output flags.Filename `short:"o" long:"output" required:"yes" value-name:"OUTPUT" description:"output destination of the alignment results" default:"-"`
}

var parser = flags.NewParser(nil, flags.Default)

func main() {
	if _, err := parser.Parse(); err != nil {
		if flagsErr, ok := err.(*flags.Error); ok && flagsErr.Type == flags.ErrHelp {
			os.Exit(0)
		} else {
			os.Exit(127)
		}
	}
}

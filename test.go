package main

import (
	"github.com/chrislusf/glow/flow"
	"math"
)

func main() {
	inChan := make(chan func() float64)
	outChan := make(chan float64)
	go func() {
		inChan <- func() float64 { return 15.1 }
		inChan <- func() float64 { return 112.5 }
		inChan <- func() float64 { return 46.3 }
		inChan <- func() float64 { return 3.7 }
		close(inChan)
	}()

	f := flow.New().
		Channel(inChan).
		Map(func(cb func() float64) float64 {
			return cb()
		}).
		Reduce(math.Max).
		AddOutput(outChan)
	go f.Run()

	for score := range outChan {
		println(score)
	}
}

package alignment

type Gene string

type PositionalIndelScores map[Gene]map[int]([2]int)

type AlignmentProfile struct {
	StopCodonPenalty         int
	GapOpenPenalty           int
	IndelCodonOpeningBonus   int
	IndelCodonExtensionBonus int
	PositionalIndelScores    map[Gene]PositionalIndelScores
}

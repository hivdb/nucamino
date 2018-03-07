package hiv1b

import (
	ap "github.com/hivdb/nucamino/alignmentprofile"
	a "github.com/hivdb/nucamino/types/amino"
)

var hiv1bPositionalIndelScores = ap.GenePositionalIndelScores{
	"GAG": map[int][2]int{
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
	"POL": map[int][2]int{
		// 56prePR + 99PR = 155
		155 + 63:  [2]int{-5, 0},
		-155 - 63: [2]int{-5, 0}, // deletion penalty to the far end of RT69
		155 + 64:  [2]int{-5, 0},
		-155 - 64: [2]int{-5, 0}, // deletion penalty to the far end of RT69
		155 + 65:  [2]int{-7, 0},
		155 + 66:  [2]int{-7, 0},
		155 + 67:  [2]int{-7, 0},
		155 + 68:  [2]int{-3, 0},
		-155 - 68: [2]int{0, 0},   // remove deletion bonus from RT68/POL223
		155 + 69:  [2]int{18, -3}, // group all insertions to RT69/POL224
		155 + 70:  [2]int{-3, 0},
		155 + 71:  [2]int{-3, 0},
		155 + 72:  [2]int{-3, 0},
		155 + 73:  [2]int{-3, 0},
	},
}

var (
	HIV1BSEQ_GAG = a.ReadString(`
		MGARASVLSGGELDRWEKIRLRPGGKKKYKLKHIVWASRELERFAVNPGLLETSEGCRQILGQLQPSLQT
		GSEELRSLYNTVATLYCVHQRIEVKDTKEALEKIEEEQNKSKKKAQQAAADTGNSSQVSQNYPIVQNLQG
		QMVHQAISPRTLNAWVKVVEEKAFSPEVIPMFSALSEGATPQDLNTMLNTVGGHQAAMQMLKETINEEAA
		EWDRLHPVHAGPIAPGQMREPRGSDIAGTTSTLQEQIGWMTNNPPIPVGEIYKRWIILGLNKIVRMYSPT
		SILDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPDCKTILKALGPAATLEEMMTAC
		QGVGGPGHKARVLAEAMSQVTNSATIMMQRGNFRNQRKTVKCFNCGKEGHIAKNCRAPRKKGCWKCGKEG
		HQMKDCTERQANFLGKIWPSHKGRPGNFLQSRPEPTAPPEESFRFGEETTTPSQKQEPIDKELYPLASLR
		SLFGNDPSSQ
	`) // 500 AAs
	HIV1BSEQ_POL = a.ReadString(`
		FFREDLAFPQGKAREFSSEQTRANSPTRRELQVWGRDNNSLSEAGADRQGTVSFSFPQITLWQRPLVTI
		KIGGQLKEALLDTGADDTVLEEMNLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVN
		IIGRNLLTQIGCTLNFPISPIETVPVKLKPGMDGPKVKQWPLTEEKIKALVEICTEMEKEGKISKIGPEN
		PYNTPVFAIKKKDSTKWRKLVDFRELNKRTQDFWEVQLGIPHPAGLKKKKSVTVLDVGDAYFSVPLDKDF
		RKYTAFTIPSINNETPGIRYQYNVLPQGWKGSPAIFQSSMTKILEPFRKQNPDIVIYQYMDDLYVGSDLE
		IGQHRTKIEELRQHLLRWGFTTPDKKHQKEPPFLWMGYELHPDKWTVQPIVLPEKDSWTVNDIQKLVGKL
		NWASQIYAGIKVKQLCKLLRGTKALTEVIPLTEEAELELAENREILKEPVHGVYYDPSKDLIAEIQKQGQ
		GQWTYQIYQEPFKNLKTGKYARMRGAHTNDVKQLTEAVQKIATESIVIWGKTPKFKLPIQKETWEAWWTE
		YWQATWIPEWEFVNTPPLVKLWYQLEKEPIVGAETFYVDGAANRETKLGKAGYVTDRGRQKVVSLTDTTN
		QKTELQAIHLALQDSGLEVNIVTDSQYALGIIQAQPDKSESELVSQIIEQLIKKEKVYLAWVPAHKGIGG
		NEQVDKLVSAGIRKVLFLDGIDKAQEEHEKYHSNWRAMASDFNLPPVVAKEIVASCDKCQLKGEAMHGQV
		DCSPGIWQLDCTHLEGKIILVAVHVASGYIEAEVIPAETGQETAYFLLKLAGRWPVKTIHTDNGSNFTST
		TVKAACWWAGIKQEFGIPYNPQSQGVVESMNKELKKIIGQVRDQAEHLKTAVQMAVFIHNFKRKGGIGGY
		SAGERIVDIIATDIQTKELQKQITKIQNFRVYYRDSRDPLWKGPAKLLWKGEGAVVIQDNSDIKVVPRRK
		AKIIRDYGKQMAGDDCVASRQDED
	`) // 56prePR + 99PR + 560RT + 288IN
	HIV1BSEQ_GP41 = a.ReadString(`
		AVGIGAMFLGFLGAAGSTMGAASMTLTVQARQLLSGIVQQQNNLLRAIEAQQHLLQLTVWGIKQLQARVL
		AVERYLKDQQLLGIWGCSGKLICTTAVPWNASWSNKSLDEIWDNMTWMEWEREIDNYTSLIYTLIEESQN
		QQEKNEQELLELDKWASLWNWFDITNWLWYIKIFIMIVGGLVGLRIVFAVLSIVNRVRQGYSPLSFQTRL
		PAPRGPDRPEGIEEEGGERDRDRSGRLVDGFLALIWDDLRSLCLFSYHRLRDLLLIVTRIVELLGRRGWE
		VLKYWWNLLQYWSQELKNSAVSLLNATAIAVAEGTDRVIEVVQRACRAILHIPRRIRQGLERALL
	`) // 345 AAs
)

var HIV1BRefLookup = ap.ReferenceSeqs{
	"GAG":  HIV1BSEQ_GAG,
	"POL":  HIV1BSEQ_POL,
	"GP41": HIV1BSEQ_GP41,
}

var Profile = ap.AlignmentProfile{
	StopCodonPenalty:         4,
	GapOpeningPenalty:        10,
	GapExtensionPenalty:      2,
	IndelCodonOpeningBonus:   0,
	IndelCodonExtensionBonus: 2,
	GeneIndelScores:          hiv1bPositionalIndelScores,
	ReferenceSequences:       HIV1BRefLookup,
}

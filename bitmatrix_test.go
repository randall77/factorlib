package factorlib

import (
	"testing"
	"sort"
)

func TestBitMatrix(t *testing.T) {
	m := newBitMatrix(20)
	
	idlist := m.addRow([]uint{1,3,5}, "A")
	if idlist != nil {
		t.Fatalf("premature result %v", idlist)
	}
	idlist = m.addRow([]uint{3,5,7}, "B")
	if idlist != nil {
		t.Fatalf("premature result %v", idlist)
	}
	idlist = m.addRow([]uint{5, 7, 9}, "C")
	if idlist != nil {
		t.Fatalf("premature result %v", idlist)
	}
	idlist = m.addRow([]uint{5, 9, 1, 11}, "D")
	if idlist != nil {
		t.Fatalf("premature result %v", idlist)
	}
	idlist = m.addRow([]uint{5, 9, 1}, "E")
	if idlist == nil {
		t.Fatalf("bad nil result")
	}
	var a []string
	for _, id := range idlist {
		a = append(a, id.(string))
	}
	sort.Strings(a)
	if a[0] != "A" || a[1] != "B" || a[2] != "C" || a[3] != "E" {
		t.Fatalf("bad ids %v", a)
	}
}

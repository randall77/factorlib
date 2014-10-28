package linear

import (
	"testing"
	"sort"
)

func TestMatrix(t *testing.T) {
	m := NewMatrix(20)
	
	idlist := m.AddRow([]uint{1,3,5}, "A")
	if idlist != nil {
		t.Fatalf("premature result %v", idlist)
	}
	idlist = m.AddRow([]uint{3,5,7}, "B")
	if idlist != nil {
		t.Fatalf("premature result %v", idlist)
	}
	idlist = m.AddRow([]uint{5, 7, 9}, "C")
	if idlist != nil {
		t.Fatalf("premature result %v", idlist)
	}
	idlist = m.AddRow([]uint{5, 9, 1, 11}, "D")
	if idlist != nil {
		t.Fatalf("premature result %v", idlist)
	}
	if m.Rows() != 4 {
		t.Fatalf("bad row count %v", idlist)
	}
	idlist = m.AddRow([]uint{5, 9, 1}, "E")
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

	m = NewMatrix(1000)
	idlist = m.AddRow([]uint{33,33,77,77}, "A")
	if idlist == nil {
		t.Fatalf("trivial row not returned")
	}
	if len(idlist) != 1 {
		t.Fatalf("trivial row length is not 1")
	}
	if idlist[0].(string) != "A" {
		t.Fatalf("trivial row did not return its id")
	}
}

package factorlib

// Implementation of a matrix over GF(2).  Used to find
// linear combinations of rows which are zero.
type bitMatrix struct {
	// size of matrix
	n int

	// ids of rows (in the order in which they were added)
	ids []interface{}

	// matrix rows
	rows []row
}

// Return a new matrix which can handle indexes 0 <= i < n.
func newBitMatrix(n int) *bitMatrix {
	return &bitMatrix{n, make([]interface{}, 0, n), make([]row, 0, n)}
}

type row struct {
	// A bit vector of 2*n bits.  The first n bits are a combination
	// of original rows.  The second n bits mark which original rows
	// were combined to make this one.
	bits bitVec

	// The column that we pivot with. == bits.firstBit()
	pivot int
}

// Adds the vector with the given set indexes (indexes may appear multiple
// times - an index is set if it appears an odd number of times).
// If there is a linear combination of the added rows that xor to the zero
// vector, addRow returns the identities of those vectors.  Otherwise returns nil.
func (m *bitMatrix) addRow(idxs []int, id interface{}) []interface{} {
	m.ids = append(m.ids, id)
	bits := newBitVec(2*m.n+1)
	for _, i := range idxs {
		bits.toggleBit(i)
	}
	bits.setBit(m.n+len(m.rows))
	for _, r := range m.rows {
		if bits.getBit(r.pivot) {
			bits.xor(r.bits)
		}
	}
	p := bits.firstBit()
	if p < m.n {
		m.rows = append(m.rows, row{bits, p})
		return nil
	}

	// found a linear combination of vectors that generates the 0 vector.
	a := []interface{}{id}
	for i := 0; i < m.n; i++ {
		if bits.getBit(m.n + i) {
			a = append(a, m.ids[i])
		}
	}
	// Note: we don't add this vector to the list of rows, as it
	// is linearly dependent.
	return a
}

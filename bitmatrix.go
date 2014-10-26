package factorlib

type bitMatrix struct {
	// size of matrix (n x n)
	n int

	// matrix rows
	rows []row
}

func newBitMatrix(n int) *bitMatrix {
	return &bitMatrix{n, make([]row, 0, n)}
}

type row struct {
	// the id associated with this row
	id interface{}

	// Which original rows were xored together to make this row
	rows bitVec

	// Xor of the original rows listed above
	bits bitVec

	// The column that we pivot with. == bits.firstBit()
	pivot int
}

// Adds the vector with the given set indexes (indexes may appear multiple
// times - an index is set if it appears an odd number of times).
// If there is a linear combination of the added rows that xor to the zero
// vector, addRow returns the identities of those vectors.  Otherwise returns nil.
func (m *bitMatrix) addRow(idxs []int, id interface{}) []interface{} {
	bits := newBitVec(m.n)
	for _, i := range idxs {
		bits.toggleBit(i)
	}
	rows := newBitVec(m.n)
	for _, r := range m.rows {
		if bits.getBit(r.pivot) {
			bits.xor(r.bits)
			rows.xor(r.rows)
		}
	}
	if bits.empty() {
		// found a linear combination of vectors that generates the 0 vector.
		a := []interface{}{id}
		for i, r := range m.rows {
			if rows.getBit(i) {
				a = append(a, r.id)
			}
		}
		// Note: we don't add this vector to the list of rows, as it
		// is linearly dependent.
		return a
	}
	rows.setBit(len(m.rows))
	m.rows = append(m.rows, row{id, rows, bits, bits.firstBit()})
	return nil
}

package linear

// Implementation of a matrix over GF(2).  Used to find
// linear combinations of rows which are zero.
type Matrix struct {
	// size of matrix
	n uint

	// ids of rows (in the order in which they were added)
	ids []interface{}

	// matrix rows
	rows []row
}

type row struct {
	// A bit vector of 2*n bits.  The first n bits are a combination
	// of original rows.  The second n bits mark which original rows
	// were combined to make this one.
	bits vector

	// The column that we pivot with. == bits.firstBit()
	pivot uint
}

// Return a new matrix which can handle indexes 0 <= i < n.
func NewMatrix(n uint) *Matrix {
	return &Matrix{n, make([]interface{}, 0, n), make([]row, 0, n)}
}

func (m *Matrix) Rows() uint {
	return uint(len(m.rows))
}

// Adds the vector with the given set indexes (indexes may appear multiple
// times - an index is set if it appears an odd number of times).
// If there is a linear combination of the added rows that xor to the zero
// vector, addRow returns the identities of those vectors.  Otherwise returns nil.
func (m *Matrix) AddRow(idxs []uint, id interface{}) []interface{} {
	bits := newVector(2*m.n)
	for _, i := range idxs {
		bits.toggleBit(i)
	}
	if bits.empty() {
		// we've been passed the all-zero vector
		return []interface{}{id}
	}
	for _, r := range m.rows {
		if bits.getBit(r.pivot) {
			bits.xor(r.bits)
		}
	}
	p := bits.firstBit()
	if p < m.n {
		bits.setBit(m.n+uint(len(m.ids)))
		m.ids = append(m.ids, id)
		m.rows = append(m.rows, row{bits, p})
		return nil
	}

	// found a linear combination of vectors that generates the 0 vector.
	a := []interface{}{id}
	for i := uint(0); i < m.n; i++ {
		if bits.getBit(m.n + i) {
			a = append(a, m.ids[i])
		}
	}
	// Note: we don't add this vector to the list of rows, as it
	// is linearly dependent.
	return a
}

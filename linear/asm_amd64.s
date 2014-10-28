#define NOSPLIT 4

// func xor(x, y *uint, size int)
// *x ^= *y, size times.
TEXT ·xor(SB),NOSPLIT,$0
	MOVQ	x+0(FP), AX
	MOVQ	y+8(FP), BX
	MOVQ	size+16(FP), CX

loop:
	CMPQ	CX, $16
	JB	tail
	// Note: we take advantage of the (undocumented) fact
	// that >= 128 byte vectors will be 16-byte aligned.
	MOVO	(AX), X0
	MOVO	16(AX), X1
	MOVO	32(AX), X2
	MOVO	48(AX), X3
	MOVO	64(AX), X4
	MOVO	80(AX), X5
	MOVO	96(AX), X6
	MOVO	112(AX), X7
	PXOR	(BX), X0
	PXOR	16(BX), X1
	PXOR	32(BX), X2
	PXOR	48(BX), X3
	PXOR	64(BX), X4
	PXOR	80(BX), X5
	PXOR	96(BX), X6
	PXOR	112(BX), X7
	MOVO	X0, (AX)
	MOVO	X1, 16(AX)
	MOVO	X2, 32(AX)
	MOVO	X3, 48(AX)
	MOVO	X4, 64(AX)
	MOVO	X5, 80(AX)
	MOVO	X6, 96(AX)
	MOVO	X7, 112(AX)
	ADDQ	$128, AX
	ADDQ	$128, BX
	SUBQ	$16, CX
	JMP	loop
tail:
	TESTQ	CX, CX
	JE	done
	MOVQ	(BX), DX
	XORQ	DX, (AX)
	ADDQ	$8, AX
	ADDQ	$8, BX
	DECQ	CX
	JMP	tail
done:
	RET

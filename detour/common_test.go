package detour

import (
	"testing"

	"github.com/stretchr/testify/assert"
)

// 중요한 알고리즘 단위 테스트를 진행하면서 내부 구현을 익힌다.
func TestOverlapBounds(t *testing.T) {
	amin := [3]float32{0, 0, 0}
	amax := [3]float32{1, 1, 1}
	bmin := [3]float32{0, 0, 0}
	bmax := [3]float32{2, 2, 2}

	r1 := OverlapBounds(amin[:], amax[:], bmin[:], bmax[:])
	assert.True(t, r1)

	amin = [3]float32{0, 0, 0}
	amax = [3]float32{1, 1, 1}
	bmin = [3]float32{1.2, 1.2, 1.2}
	bmax = [3]float32{2, 2, 2}

	r1 = OverlapBounds(amin[:], amax[:], bmin[:], bmax[:])
	assert.True(t, !r1)
}

func TestILog2(t *testing.T) {

}

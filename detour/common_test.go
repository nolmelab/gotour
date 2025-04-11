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

// 밑 2의 정수 로그
func TestILog2(t *testing.T) {
	assert.True(t, Ilog2(3) == 1)
	assert.True(t, Ilog2(127) == 6)
}

// 벡터 함수들
func TestVec(t *testing.T) {
	v1 := []float32{0, 0, 0}
	v2 := []float32{1, 0, 0}
	assert.True(t, Vdist(v1, v2) == 1)
}

func TestClosestPointTriangle(t *testing.T) {
	// AI시대이다. Grok은 코드를 단계별로 상세하게 설명한다. Gemini는 수학공식을 보여준다.
	a := []float32{0, 0, 0}
	b := []float32{0, 1, 0}
	c := []float32{1, 1, 1}
	p := []float32{10, 10, 10}

	var closest [3]float32

	ClosestPtPointTriangle(closest[:], p, a, b, c)
	t.Log(closest)

	// 계산 기하의 테스트는 쉽지 않은 작업이다. 정확하게 계산하고 그 값들을 사용해야 한다.
	// case들을 늘려 가도록 한다.
}

func TestClosestHeightPointTriange(t *testing.T) {

}
